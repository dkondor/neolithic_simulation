/*  -*- C++ -*-
 * read_table_cpp.h -- simple and robust general methods for reading numeric data
 * 	from text files, e.g. TSV or CSV
 * 
 * simple: should be usable in a few lines of code
 * robust: try to detect and signal errors (in format, overflow, underflow etc.)
 * 	especially considering cases that would be silently ignored with
 * 	scanf / iostreams or similar; avoid undefined behavior
 * 
 * C++-only version, which does not require the POSIX getline() function,
 * uses std::getline() from the C++ standard library which should be available
 * on all platforms
 * 
 * note that this requires C++11
 * 
 * Copyright 2018 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 * example usage in C++

istream f; ... // open input file or use stdin
read_table2 r(f);
while(r.read_line()) {
	int32_t id1,id2;
	double d1;
	uint64_t id3;
	if(!r.read(id1,d1,id3,id2)) break; // false return value indicates error
	... // do something with the values read
}
if(r.get_last_error() != T_EOF) { // handle error
	std::cerr<<"Error reading input:\n";
	r.write_error(std::cerr);
}

 */

#ifndef _READ_TABLE_H
#define _READ_TABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>
#include <utility>
#include <type_traits>
#include <memory>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#if __cplusplus >= 201703L
#include <string_view>
#endif

#include <cmath>

/* possible error codes */
enum read_table_errors {T_OK = 0, T_EOF = 1, T_EOL, T_MISSING, T_FORMAT,
	T_OVERFLOW, T_NAN, T_TYPE, T_COPIED, T_ERROR_FOPEN, T_READ_ERROR};
static const char * const error_desc[] = {"No error", "End of file", "Unexpected end of line",
		"Missing value", "Invalid value", "Overflow or underflow", "NaN or infinity read",
		"Unknown conversion requested", "Invalidated instance", "Error opening file",
		"Error reading input"};

/* convert error code to string description */
static const char* get_error_desc(enum read_table_errors err) {
	const static char unkn[] = "Unknown error";
	switch(err) {
		case T_OK:
			return error_desc[0];
		case T_EOF:
			return error_desc[1];
		case T_EOL:
			return error_desc[2];
		case T_MISSING:
			return error_desc[3];
		case T_FORMAT:
			return error_desc[4];
		case T_OVERFLOW:
			return error_desc[5];
		case T_NAN:
			return error_desc[6];
		case T_TYPE:
			return error_desc[7];
		case T_COPIED:
			return error_desc[8];
		case T_ERROR_FOPEN:
			return error_desc[9];
		case T_READ_ERROR:
			return error_desc[10];
		default:
			return unkn;
	}
}


/* helper classes to be given as parameters to read_table2::read_next() and
 * read_table2::read() */
 
 
/* dummy struct to be able to call the same interface to skip data */
struct read_table_skip_t { };
static const read_table_skip_t _read_table_skip1;
inline const read_table_skip_t& read_table_skip() { return _read_table_skip1; }
/* helper class to return read-only string parts
 * (if std::string_view is not available) */
struct string_view_custom {
	const char* str;
	size_t len;
	const char* data() const { return str; }
	size_t length() const { return len; }
	size_t size() const { return len; }
	char operator [] (size_t i) const { return str[i]; }
	int print(FILE* f) const {
		if(len == 0) return 0;
		if(len <= INT32_MAX) return fprintf(f,"%.*s",(int)len,str);
		return -1;
	}
	string_view_custom():str(0),len(0) { }
	bool operator == (const string_view_custom& v) const {
		if(len != v.len) return false; /* lengths must be the same */
		if(len == 0) return true; /* empty strings are considered equal */
		if(str && v.str) return strncmp(str,v.str,len) == 0;
		else return false; /* str or v.str is null, this is probably an error */
	}
};
template<class ostream>
ostream& operator << (ostream& s, const string_view_custom& str) {
	s.write(str.str,str.len);
	return s;
}

/* struct to represent values with bounds */
template<class T>
struct read_bounds_t {
	read_bounds_t(T& val_, T min_, T max_):val(val_),min(min_),max(max_) { }
	T& val;
	T min;
	T max;
};
template<class T> read_bounds_t<T> read_bounds(T& val_, T min_, T max_) {
	return read_bounds_t<T>(val_,min_,max_);
}
/* shortcut to read coordinate pairs in the "obvious" format, i.e. the first
 * value should be between -180.0 and 180.0, the second value should be
 * between -90.0 and 90.0
 * -- note: this is the format that is obvious to me, different use cases
 * might have the coordinates in different order or the range of longitudes
 * could be 0 to 360 or even unbounded -- feel free to modify the bounds :) */
inline read_bounds_t<std::pair<double,double> > read_bounds_coords(std::pair<double,double>& coords) {
	return read_bounds_t<std::pair<double,double> >(coords,
		std::make_pair(-180.0,-90.0),std::make_pair(180.0,90.0));
}

struct line_parser_params {
	int base; /* base for integer conversions */
	char delim; /* delimiter to use; 0 means any blank (space or tab) note: cannot be newline */
	char comment; /* character to indicate comments; 0 means none */
	bool allow_nan_inf; /* further flags: whether reading a NaN or INF for double values is considered and error */
	line_parser_params():base(10),delim(0),comment(0),allow_nan_inf(true) { }
	line_parser_params& set_base(int base_) { base = base_; return *this; }
	line_parser_params& set_delim(char delim_) { delim = delim_; return *this; }
	line_parser_params& set_comment(char comment_) { comment = comment_; return *this; }
	line_parser_params& set_allow_nan_inf(bool allow_nan_inf_) { allow_nan_inf = allow_nan_inf_; return *this; }
};

/* "helper" class doing most of the work for parsing only one line */
struct line_parser {
	protected:
		std::string buf; /* buffer to hold the current line */
		size_t pos = 0; /* current position in line */
		size_t col = 0; /* current field (column) */
		int base; /* base for integer conversions */
		enum read_table_errors last_error = T_OK; /* error code of the last operation */
		char delim; /* delimiter to use; 0 means any blank (space or tab) note: cannot be newline */
		char comment; /* character to indicate comments; 0 means none */
		bool allow_nan_inf; /* further flags: whether reading a NaN or INF for double values is considered and error */
		
		void line_parser_init(const line_parser_params& par) {
			base = par.base;
			delim = par.delim;
			comment = par.comment;
			allow_nan_inf = par.allow_nan_inf;
		}
		
	public:
		/* 1. constructors, either with an empty string, or anything that can be copied into a string */
		explicit line_parser(const line_parser_params& par = line_parser_params()) {
			line_parser_init(par);
		}
		template<class... Args>
		explicit line_parser(const line_parser_params& par, Args&&... args):buf(std::forward<Args>(args)...) {
			line_parser_init(par);
		}
		template<class T, class... Args, std::enable_if_t<!std::is_base_of<line_parser, T>::value>* = nullptr>
		explicit line_parser(T&& t, Args&&... args):buf(std::forward<T>(t), std::forward<Args>(args)...) {
			line_parser_init(line_parser_params());
		}
		
		line_parser(const line_parser& lp) = default;
		line_parser& operator = (const line_parser& lp) = default;
		
		/* move constructor and move assignment -- ensure the string is moved */
		line_parser(line_parser&& lp) : buf(std::move(lp.buf)) {
			pos = lp.pos;
			col = lp.col;
			base = lp.base;
			comment = lp.comment;
			allow_nan_inf = lp.allow_nan_inf;
			last_error = lp.last_error;
			lp.last_error = T_COPIED;
			lp.pos = 0;
			lp.col = 0;
		}
		line_parser& operator = (line_parser&& lp) {
			if(this == &lp) return *this; /* protect self-assignment */
			buf = std::move(lp.buf);
			pos = lp.pos;
			col = lp.col;
			base = lp.base;
			comment = lp.comment;
			allow_nan_inf = lp.allow_nan_inf;
			last_error = lp.last_error;
			lp.last_error = T_COPIED;
			lp.pos = 0;
			lp.col = 0;
			return *this;
		}
		
		/* 2. set (copy) the internal string */
		template<class... Args>
		void set_line(Args&&... args) {
			buf = std::string(std::forward<Args>(args)...);
			col = 0;
			pos = 0;
			last_error = T_OK;
		}
		/* get current line string */
		const char* get_line_c_str() const { return buf.c_str(); }
		const std::string& get_line_str() const { return buf; }
		
		/* 3. main interface for parsing data; this uses templates and has 
		 * specializations for all data types supported:
		 * 16, 32 and 64-bit signed and unsigned integers, doubles,
		 * strings, pairs of doubles and special "types":
		 * 	- read_table_skip_t for skipping values
		 * 	- read_bounds_t for specifying minimum and maximum value for the input
		 * see below for more explanation */
		/* try to parse one value from the currently read line */
		template<class T> bool read_next(T& val, bool advance_pos = true);
		/* overload of the previous for reading values with bounds */
		template<class T> bool read_next(read_bounds_t<T> val, bool advance_pos = true);
		/* try to parse whole line (read previously with read_line()),
		 * into the given list of parameters */
		bool read() { return true; }
		template<class first, class ...rest>
		bool read(first&& val, rest&&... vals);
		
		
		/* 4. functions for setting parameters */
		/* set delimiter character (default is spaces and tabs) */
		void set_delim(char delim_) { delim = delim_; }
		/* get delimiter character (default is spaces and tabs) */
		char get_delim() const { return delim; }
		/* set comment character (default is none) */
		void set_comment(char comment_) { comment = comment_; }
		/* get comment character (default is none) */
		char get_comment() const { return comment; }
		line_parser_params get_params() const {
			return line_parser_params().set_base(base).set_delim(delim).set_allow_nan_inf(allow_nan_inf).set_comment(comment);
		}
		void reset_pos() {
			if(last_error == T_COPIED || last_error == T_EOF ||
				last_error == T_ERROR_FOPEN || last_error == T_READ_ERROR)
					return; /* these errors cannot be recovered, and should not be reset */
			
			pos = 0;
			col = 0;
			last_error = T_OK;
		}
		
		/* get last error code */
		enum read_table_errors get_last_error() const { return last_error; }
		const char* get_last_error_str() const { return get_error_desc(last_error); }
		
		size_t get_pos() const { return pos; }
		size_t get_col() const { return col; }
		
		
		/* 5. non-templated functions for reading specific data types and values */
		/* skip next field, ignoring any content */
		bool read_skip();
		/* read one 32-bit signed integer in the given limits */
		bool read_int32_limits(int32_t& i, int32_t min, int32_t max, bool advance_pos = true);
		bool read_int32(int32_t& i, bool advance_pos = true) { return read_int32_limits(i,INT32_MIN,INT32_MAX,advance_pos); }
		/* read one 32-bit unsigned integer in the given limits */
		bool read_uint32_limits(uint32_t& i, uint32_t min, uint32_t max, bool advance_pos = true);
		bool read_uint32(uint32_t& i, bool advance_pos = true) { return read_uint32_limits(i,0,UINT32_MAX,advance_pos); }
		/* read one 64-bit signed integer in the given limits */
		bool read_int64_limits(int64_t& i, int64_t min, int64_t max, bool advance_pos = true);
		bool read_int64(int64_t& i, bool advance_pos = true) { return read_int64_limits(i,INT64_MIN,INT64_MAX,advance_pos); }
		/* read one 64-bit unsigned integer in the given limits */
		bool read_uint64_limits(uint64_t& i, uint64_t min, uint64_t max, bool advance_pos = true);
		bool read_uint64(uint64_t& i, bool advance_pos = true) { return read_uint64_limits(i,0,UINT64_MAX,advance_pos); }
		/* read one 16-bit signed integer in the given limits */
		bool read_int16_limits(int16_t& i, int16_t min, int16_t max, bool advance_pos = true);
		bool read_int16(int16_t& i, bool advance_pos = true) { return read_int16_limits(i,INT16_MIN,INT16_MAX,advance_pos); }
		/* read one 16-bit unsigned integer in the given limits */
		bool read_uint16_limits(uint16_t& i, uint16_t min, uint16_t max, bool advance_pos = true);
		bool read_uint16(uint16_t& i, bool advance_pos = true) { return read_uint16_limits(i,0,UINT16_MAX,advance_pos); }
		/* read one double value in the given limits */
		bool read_double_limits(double& d, double min, double max, bool advance_pos = true);
		bool read_double(double& d, bool advance_pos = true);
		/* read string, copying from the buffer */
		bool read_string(std::string& str, bool advance_pos = true);
		/* read string, return readonly view */
#if __cplusplus >= 201703L
		bool read_string_view(std::string_view& str, bool advance_pos = true);
#endif
		bool read_string_view_custom(string_view_custom& str, bool advance_pos = true);
	
	protected:
		/* helper functions for the previous */
		bool read_table_pre_check(bool advance_pos);
		bool read_table_post_check(const char* c2);
		/* read string return start position and length
		 *  -- the other read_string functions then use these to create the string_view or copy to a string */
		bool read_string2(std::pair<size_t,size_t>& pos1, bool advance_pos = true);
};


/* main class containing main parameters for processing text */
struct read_table2 : public line_parser {
	protected:
		std::istream* is = nullptr; /* input stream -- note: only a pointer is stored, the caller either supplies an input stream or a
			file name; in the former case, the original object should not go out of scope while this struct is used */
		std::unique_ptr<std::ifstream> fs; /* file stream if it is opened by us */
		const char* fn = nullptr; /* file name, stored optionally for error output (note: not owned by this class, caller should not free the supplied value) */
		uint64_t line = 0; /* current line (count starts from 1) */
		read_table2() = delete; /* user should supply either an input stream or a filename to use */
		read_table2(const read_table2& r) = delete; /* disallow copying, only moving is possible */
		
		const char line_endings[2] = {'\n','\r'};
		
	public:
		
		/* 1. constructors -- need to give a file name or an already open input stream */
		
		/* constructor taking a file name; the given file is opened for reading
		 * and closed later in the destructor */
		explicit read_table2(const char* fn_, const line_parser_params& par = line_parser_params());
		/* constructor taking a reference to a stream -- it is NOT copied, i.e.
		 * the original instance of the stream need to be kept by the caller;
		 * also, the stream is not closed in the destructor in this case */
		explicit read_table2(std::istream& is_, const line_parser_params& par = line_parser_params());
		/* constructor either opening a file or taking the stream as fallback
		 * when fn_ is NULL */
		read_table2(const char* fn_, std::istream& is_, const line_parser_params& par = line_parser_params());
		read_table2(read_table2&& r);
		
		read_table2& operator = (read_table2&& r);
		
		/* 2. read one line into the internal buffer
		 * 	the 'skip' parameter controls whether empty lines are skipped */
		bool read_line(bool skip = true);
		
		/* get current position in the file */
		uint64_t get_line() const { return line; }
		/* set filename (for better formatting of diagnostic messages) */
		void set_fn_for_diag(const char* fn_) { fn = fn_; }
		const char* get_fn() const { return fn; }
		
		/* write formatted error message to the given stream */
		void write_error(std::ostream& f) const;
		void write_error(FILE* f) const;
		
		/* create a string error message that can be thrown as an exception */
		std::string exception_string(std::string&& base_message = "") {
			std::ostringstream strs(std::move(base_message), std::ios_base::ate);
			strs << "read_table, ";
			if(fn) strs << "file " << fn << ", ";
			else strs << "input ";
			strs << "line " << line << ", position " << pos << " / column " << col << ": ";
			strs << get_error_desc(last_error) << '\n';
			return strs.str();
		}
};


/* templated functions, template specializations to automatically handle
 * reading the proper types */

/* template specializations to use the same function name */
template<> inline bool line_parser::read_next(int32_t& val, bool advance_pos) { return read_int32(val,advance_pos); }
template<> inline bool line_parser::read_next(uint32_t& val, bool advance_pos) { return read_uint32(val,advance_pos); }
template<> inline bool line_parser::read_next(int16_t& val, bool advance_pos) { return read_int16(val,advance_pos); }
template<> inline bool line_parser::read_next(uint16_t& val, bool advance_pos) { return read_uint16(val,advance_pos); }
template<> inline bool line_parser::read_next(int64_t& val, bool advance_pos) { return read_int64(val,advance_pos); }
template<> inline bool line_parser::read_next(uint64_t& val, bool advance_pos) { return read_uint64(val,advance_pos); }
template<> inline bool line_parser::read_next(double& val, bool advance_pos) { return read_double(val,advance_pos); }
template<> inline bool line_parser::read_next(std::pair<double,double>& p, bool advance_pos) {
	double x,y;
	size_t old_pos = pos;
	bool ret = read_double(x) && read_double(y);
	if(ret) p = std::pair<double,double>(x,y);
	if(!advance_pos) pos = old_pos;
	return ret;
}
template<> inline bool line_parser::read_next(std::string& str, bool advance_pos) { return read_string(str,advance_pos); }
#if __cplusplus >= 201703L
template<> inline bool line_parser::read_next(std::string_view& str, bool advance_pos) { return read_string_view(str,advance_pos); }
#endif
template<> inline bool line_parser::read_next(string_view_custom& str, bool advance_pos) { return read_string_view_custom(str,advance_pos); }

/* dummy struct to be able to call the same interface to skip data
 * (useful if used with the variadic template below) */
template<> inline bool line_parser::read_next(const read_table_skip_t&, bool) { return read_skip(); }
//~ template<> bool line_parser::read_next(read_table_skip_t skip) { return read_skip(); }


/* overloads for reading with bounds
 * example usage:
line_parser r(...);
uint32_t x;
r.read_next(read_bounds(x,1000U,2000U));
*/
template<> inline bool line_parser::read_next(read_bounds_t<int32_t> b, bool advance_pos) {
	return read_int32_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<uint32_t> b, bool advance_pos) {
	return read_uint32_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<int64_t> b, bool advance_pos) {
	return read_int64_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<uint64_t> b, bool advance_pos) {
	return read_uint64_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<int16_t> b, bool advance_pos) {
	return read_int16_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<uint16_t> b, bool advance_pos) {
	return read_uint16_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<double> b, bool advance_pos) {
	return read_double_limits(b.val,b.min,b.max,advance_pos);
}
template<> inline bool line_parser::read_next(read_bounds_t<std::pair<double,double> > b, bool advance_pos) {
	double x,y;
	size_t old_pos = pos;
	bool ret = read_double_limits(x,b.min.first,b.max.first) &&
		read_double_limits(y,b.min.second,b.max.second);
	if(ret) b.val = std::make_pair(x,y);
	if(!advance_pos) pos = old_pos;
	return ret;
}

/* recursive templated function to convert whole line using one function call only
 * note: recursion will be probably eliminated and the whole function expanded to
 * the actual sequence of conversions needed */
//~ bool line_parser::read() { return true; }
template<class first, class ...rest>
bool line_parser::read(first&& val, rest&&... vals) {
	if(!read_next(val,true)) return false;
	return read(vals...);
}


#endif /* _READ_TABLE_H */

