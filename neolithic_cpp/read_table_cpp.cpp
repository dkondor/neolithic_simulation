/*
 * read_table_cpp.cpp -- simple and robust general methods for reading numeric data
 * 	from text files, e.g. TSV or CSV
 * 
 * C++ implementation separated from the header
 * 
 * Copyright 2021 Daniel Kondor <kondor.dani@gmail.com>
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
 */


#include "read_table_cpp.h"



read_table2::read_table2(const char* fn_, const line_parser_params& par) : line_parser(par) {
	fs = std::unique_ptr<std::ifstream>(new std::ifstream(fn_));
	if( !fs || !(fs->is_open()) || fs->fail() ) last_error = T_ERROR_FOPEN;
	else fs->exceptions(std::ios_base::goodbit); /* clear exception mask -- no exceptions thrown, error checking done separately */
	is = fs.get();
	fn = fn_;
}

read_table2::read_table2(const char* fn_, std::istream& is_, const line_parser_params& par) : line_parser(par) {
	line_parser_init(par);
	if(fn_) {
		fs = std::unique_ptr<std::ifstream>(new std::ifstream(fn_));
		if( !fs || !(fs->is_open()) || fs->fail() ) last_error = T_ERROR_FOPEN;
		else fs->exceptions(std::ios_base::goodbit); /* clear exception mask -- no exceptions thrown, error checking done separately */
		is = fs.get();
	}
	else {
		is = &is_;
		is->exceptions(std::ios_base::goodbit); /* clear exception mask -- no exceptions thrown, error checking done separately */
	}
	fn = fn_;
}

read_table2::read_table2(std::istream& is_, const line_parser_params& par) : line_parser(par), is(&is_) {
	is->exceptions(std::ios_base::goodbit); /* clear exception mask -- no exceptions thrown, error checking done separately */
}


/* move constructor -- moves the stream to the new instance
 * the old instance is invalidated */
read_table2::read_table2(read_table2&& r) : line_parser(std::move(r)), 
		is(r.is), fs(std::move(r.fs)), fn(r.fn), line(r.line) {
	/* note: line_parser base class' move constructor will set r.last_error == T_COPIED,
	 * so r will not be usable from this point on */
	r.is = nullptr;
}

read_table2& read_table2::operator = (read_table2&& r) {
	if(this == &r) return *this;
	line_parser::operator =(std::move(r)); /* move the main part, also setting r.last_error == T_COPIED */ 
	is = r.is;
	fs = std::move(r.fs);
	fn = r.fn;
	line = r.line;
	r.is = nullptr;
	return *this;
}

/* read a new line (discarding any remaining data in the current line)
 * returns true if a line was read, false on failure
 * note that failure can mean end of file, which should be checked separately
 * if skip == 1, empty lines are skipped (i.e. reading continues until a
 * nonempty line is found); otherwise, empty lines are read and stored as well,
 * which will probably result in errors if data is tried to be parsed from it */
bool read_table2::read_line(bool skip) {
	if(last_error == T_EOF || last_error == T_COPIED ||
		last_error == T_ERROR_FOPEN) return false;
	if(is->eof()) { last_error = T_EOF; return false; }
	while(1) {
		std::getline(*is,buf);
		if(is->eof()) { last_error = T_EOF; return false; }
		if(is->fail()) { last_error = T_READ_ERROR; return false; }
		size_t len = buf.size();
		line++; 
		pos = 0;
		/* check that there is actual data in the line, empty lines are skipped */
		if(skip) {
			for(; pos < len; pos++)
				if( ! (buf[pos] == ' ' || buf[pos] == '\t') ) break;
			if(pos < len) {
				if(comment && buf[pos] == comment) continue; /* if there is only a comment, then this line is skipped */
				if(delim) pos = 0; /* if there is a delimiter character then whitespace at the beginning of a line is not skipped */
				break; /* there is some data in the line */
			}
		}
		else break; /* if empty lines should not be skipped */
	}
	col = 0; /* reset the counter for columns */
	last_error = T_OK;
	
	/* check line end characters */
	if(buf.size()) for(char c : line_endings) if(buf.back() == c) { buf.pop_back(); break; }
	
	return true;
}

/* checks to be performed before trying to convert a field */
bool line_parser::read_table_pre_check(bool advance_pos) {
	if(last_error == T_EOF || last_error == T_EOL || last_error == T_COPIED ||
		last_error == T_READ_ERROR || last_error == T_ERROR_FOPEN) return false;
	/* 1. skip any blanks */
	size_t old_pos = pos;
	size_t len = buf.size();
	for(;pos < len; pos++)
		if( ! (buf[pos] == ' ' || buf[pos] == '\t') ) break;
	/* 2. check for end of line or comment */
	if(pos == len || buf[pos] == '\n' || (comment && buf[pos] == comment) ) {
		last_error = T_EOL;
		if(!advance_pos) pos = old_pos;
		return false;
	}
	/* 3. check for field delimiter (if we have any) */
	if(delim && buf[pos] == delim) {
		last_error = T_MISSING;
		if(!advance_pos) pos = old_pos;
		return false;
	}
	return true;
}

/* perform checks needed after number conversion */
bool line_parser::read_table_post_check(const char* c2) {
	/* 0. check for format errors and overflow as indicated by strto* */
	if(errno == EINVAL || c2 == buf.c_str() + pos) {
		last_error = T_FORMAT;
		return false;
	}
	if(errno == ERANGE) {
		last_error = T_OVERFLOW;
		return false;
	}
	/* 1. skip the converted number and any blanks */
	bool have_blank = false;
	size_t len = buf.size();
	for(pos = c2 - buf.c_str();pos<len;pos++)
		if( ! (buf[pos] == ' ' || buf[pos] == '\t') ) break;
		else have_blank = true;
	last_error = T_OK;
	/* 2. check for end of line -- this is not a problem here */
	if(pos == len || buf[pos] == '\n' ||
		(comment && buf[pos] == comment) ) return true;
	if(delim == 0 && have_blank == false) {
		/* if there is no explicit delimiter, then there need to be at least
		 * one blank after the converted number if it is not the end of line */
		last_error = T_FORMAT;
		return false;
	}
	/* 3. otherwise, check for proper delimiter, if needed */
	if(delim) {
		if(buf[pos] != delim) {
			last_error = T_FORMAT;
			return false;
		}
		pos++; /* in this case, advance position further, past the delimiter */
	}
	col++; /* advance column counter as well */
	return true; /* everything OK */
}


/* skip next field, ignoring any content
 * if we have a delimiter, this means advancing until the next delimiter and
 * 	then one more position
 * if no delimiter, this means skipping any blanks, than any nonblanks and
 * 	ending at the next blank */
bool line_parser::read_skip() {
	size_t len = buf.size();
	if(delim) {
		/* if there is a delimiter, just advance until after the next one */
		for(;pos<len;pos++) if(buf[pos] == delim) break;
		if(pos == len) {
			last_error = T_EOL;
			return false;
		}
		pos++; /* note: we do not care what is after the delimiter */
	}
	else {
		/* no delimiter, skip any blanks, then skip all non-blanks */
		for(;pos<len;pos++)
			if( ! (buf[pos] == ' ' || buf[pos] == '\t') ) break;
		if(pos == len || buf[pos] == '\n' || (comment && buf[pos] == comment)) {
			last_error = T_EOL;
			return false;
		}
		for(;pos<len;pos++) if(buf[pos] == ' ' || buf[pos] == '\t' ||
			buf[pos] == '\n' || (comment && buf[pos] == comment)) break;
		/* we do not care what is after the field, now we are either at a
		* 	blank or line end */
	}
	
	col++;
	last_error = T_OK;
	return true;
}


/* return the string value in the next field -- internal helper */
bool line_parser::read_string2(std::pair<size_t,size_t>& pos1, bool advance_pos) {
	size_t len = buf.size();
	size_t old_pos = pos;
	if(delim) {
		if(last_error == T_EOF || last_error == T_EOL ||
			last_error == T_READ_ERROR || last_error == T_ERROR_FOPEN) return false;
		/* note: having an empty string is OK in this case */
		size_t p1 = pos; /* start of the string */
		for(;pos<len;pos++) if(buf[pos] == delim || buf[pos] == '\n' ||
			(comment && buf[pos] == comment)) break;
		pos1.first = p1;
		pos1.second = pos - p1;
		if(pos<len && buf[pos] == delim) pos++; /* note: we do not care what is after the delimiter */
		else last_error = T_EOL; /* save that we were already at the end of a line; trying to read another field will result in an error */
	}
	else {
		if(!read_table_pre_check(advance_pos)) return false;
		size_t p1 = pos; /* start of the string */
		for(;pos<len;pos++) if(buf[pos] == ' ' || buf[pos] == '\t' || buf[pos] == '\n' ||
			(comment && buf[pos] == comment)) break;
		/* we do not care what is after the field, now we are either at a
		 * 	blank or line end */
		pos1.first = p1;
		pos1.second = pos - p1;
	}
	if(!advance_pos) pos = old_pos;
	return true;
}

#if __cplusplus >= 201703L
/* return the string value in the next field as a string_view
 * NOTE: it will be invalidated when a new line is read */
bool line_parser::read_string_view(std::string_view& str, bool advance_pos) {
	std::pair<size_t,size_t> pos1;
	if(!read_string2(pos1,advance_pos)) return false;
	str = std::string_view(buf.data() + pos1.first, pos1.second);
	return true;
}
#endif
/* same but using a custom class instead of relying on the C++17
 * std::string_view */
bool line_parser::read_string_view_custom(string_view_custom& str, bool advance_pos) {
	std::pair<size_t,size_t> pos1;
	if(!read_string2(pos1,advance_pos)) return false;
	str.str = buf.data() + pos1.first;
	str.len = pos1.second;
	return true;
}

/* return the string value in the next field as a copy */
bool line_parser::read_string(std::string& str, bool advance_pos) {
	std::pair<size_t,size_t> pos1;
	if(!read_string2(pos1,advance_pos)) return false;
	str.assign(buf,pos1.first,pos1.second);
	return true;
}


/* try to convert the next value to integer
 * check explicitely that it is within the limits provided
 * (note: the limits are inclusive, so either min or max is OK)
 * return true on success, false on error */
bool line_parser::read_int32_limits(int32_t& i, int32_t min, int32_t max, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	bool ret;
	char* c2;
	long res = strtol(buf.c_str() + pos, &c2, base);
	/* check for format errors first (this also advances pos) */
	ret = read_table_post_check(c2);
	if(ret) {
		/* format is OK, check that result fits in bounds */
		if(res > (long)max || res < (long)min) {
			last_error = T_OVERFLOW;
			if(res > (long)max) i = max;
			if(res < (long)min) i = min;
			ret = false;
		}
		else i = res; /* store potential result */
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}

/* try to convert the next value to 64-bit integer
 * return true on success, false on error */
bool line_parser::read_int64_limits(int64_t& i, int64_t min, int64_t max, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	bool ret = true;
	char* c2;
	/* note: try to determine if we should use long or long long */
	long res;
	long long res2;
	constexpr bool use_long = (LONG_MAX >= INT64_MAX && LONG_MIN <= INT64_MIN);
	if(use_long) res = strtol(buf.c_str() + pos, &c2, base);
	else res2 = strtoll(buf.c_str() + pos, &c2, base);
	/* check for format errors and advance the position */
	ret = read_table_post_check(c2);
	if(ret) {
		/* format is OK, check for overflows */
		if(use_long) {
			if(res > (long)max || res < (long)min) {
				last_error = T_OVERFLOW;
				if(res > (long)max) i = max;
				if(res < (long)min) i = min;
				ret = false;
			}
			else i = (int64_t)res; /* store the result */
		}
		else {
			if(res2 > (long long)max || res2 < (long long)min) {
				last_error = T_OVERFLOW;
				if(res2 > (long long)max) i = max;
				if(res2 < (long long)min) i = min;
				ret = false;
			}
			else i = (int64_t)res2; /* store the result */
		}
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}

/* try to convert the next value to 32-bit unsigned integer
 * return true on success, false on error */
bool line_parser::read_uint32_limits(uint32_t& i, uint32_t min, uint32_t max, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	bool ret;
	char* c2;
	/* stricly require that the next character is alphanumeric
	 * -- strtoul() will silently accept and negate negative values */
	if( ! (isalnum(buf[pos]) || buf[pos] == '+') ) {
		if(buf[pos] == '-') last_error = T_OVERFLOW;
		else last_error = T_FORMAT;
		i = 0;
		ret = false;
	}
	else {
		unsigned long res = strtoul(buf.c_str() + pos, &c2, base);
		/* check for format errors first (this also advances pos) */
		ret = read_table_post_check(c2);
		if(ret) {
			/* format is OK, check that result fits in bounds */
			if(res > (unsigned long)max || res < (unsigned long)min) {
				last_error = T_OVERFLOW;
				if(res > (unsigned long)max) i = max;
				if(res < (unsigned long)min) i = min;
				ret = false;
			}
			else i = res; /* store potential result */
		}
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}

/* try to convert the next value to 64-bit unsigned integer
 * return true on success, false on error */
bool line_parser::read_uint64_limits(uint64_t& i, uint64_t min, uint64_t max, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	bool ret = true;
	char* c2;
	/* stricly require that the next character is alphanumeric
	 * -- strtoul() will silently accept and negate negative values */
	if( ! (isalnum(buf[pos]) || buf[pos] == '+') ) {
		if(buf[pos] == '-') last_error = T_OVERFLOW;
		else last_error = T_FORMAT;
		i = 0;
		ret = false;
	}
	else {
		/* note: try to determine if to use long or long long */
		unsigned long res;
		unsigned long long res2;
		constexpr bool use_ulong = (ULONG_MAX >= UINT64_MAX);
		if(use_ulong) res = strtoul(buf.c_str() + pos, &c2, base);
		else res2 = strtoull(buf.c_str() + pos, &c2, base);
		/* check for parse errors and advance the position */
		ret = read_table_post_check(c2);
		if(ret) {
			/* format is OK, check for overflow */
			if(use_ulong) {
				if(res > (unsigned long)max || res < (unsigned long)min) {
					last_error = T_OVERFLOW;
					if(res > (unsigned long)max) i = max;
					if(res < (unsigned long)min) i = min;
					ret = false;
				}
				else i = res; /* store the result */
			}
			else {
				if(res2 > (unsigned long long)max || res2 < (unsigned long long)min) {
					last_error = T_OVERFLOW;
					if(res2 > (unsigned long long)max) i = max;
					if(res2 < (unsigned long long)min) i = min;
					ret = false;
				}
				else i = res2; /* store the result */
			}
		}
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}

/* try to convert the next value to a 16-bit signed integer
 * return true on success, false on error
 * note: this uses the previous functions as there is no separate library
 * function for 16-bit integers anyway */
bool line_parser::read_int16_limits(int16_t& i, int16_t min, int16_t max, bool advance_pos) {
	/* just use the previous function and check for overflow */
	int32_t i2;
	/* note: the following function already check for overflow as well */
	bool ret = read_int32_limits(i2,(int32_t)min,(int32_t)max,advance_pos);
	if(ret) i = i2;
	return ret;
}

/* try to convert the next value to a 16-bit unsigned integer
 * return true on success, false on error */
bool line_parser::read_uint16_limits(uint16_t& i, uint16_t min, uint16_t max, bool advance_pos) {
	/* just use the previous function and check for overflow */
	uint32_t i2;
	bool ret = read_uint32_limits(i2,(uint32_t)min,(uint32_t)max,advance_pos);
	if(ret) i = i2;
	return ret;
}

/* try to convert the next value to a double precision float value
 * return true on success, false on error */
bool line_parser::read_double(double& d, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	char* c2;
	d = strtod(buf.c_str() + pos, &c2);
	/* advance position after the number, check if there is proper field separator */
	bool ret = read_table_post_check(c2);
	if(ret && allow_nan_inf == false) {
		if(std::isnan(d) || std::isinf(d)) {
			last_error = T_NAN;
			ret = false;
		}
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}
bool line_parser::read_double_limits(double& d, double min, double max, bool advance_pos) {
	size_t old_pos = pos;
	if(!read_table_pre_check(advance_pos)) return false;
	errno = 0;
	char* c2;
	d = strtod(buf.c_str() + pos, &c2);
	bool ret = read_table_post_check(c2);
	if(ret) {
		if(std::isnan(d)) {
			last_error = T_NAN;
			ret = false;
		}
		else {
			/* note: this will not be true if min or max is NaN, not sure if
			 * that's a problem */
			if( ! (d <= max && d >= min) ) {
				last_error = T_OVERFLOW;
				ret = false;
			}
		}
	}
	if(!advance_pos) pos = old_pos;
	return ret;
}



/* write formatted error message to the given stream */
void read_table2::write_error(std::ostream& f) const {
	f<<"read_table, ";
	if(fn) f<<"file "<<fn<<", ";
	else f<<"input ";
	f<<"line "<<line<<", position "<<pos<<" / column "<<col<<": "<<get_error_desc(last_error)<<"\n";
}

void read_table2::write_error(FILE* f) const {
	if(!f) return;
	fprintf(f,"read_table, ");
	if(fn) fprintf(f,"file %s, ",fn);
	else fprintf(f,"input ");
	fprintf(f,"line %lu, position %lu / column %lu: %s\n",line,pos,col,get_error_desc(last_error));
}



