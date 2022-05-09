/*
 * cell_it_helper.hpp -- simple iterator proxy class used by the cell
 * 	colletion classes
 * 
 * Copyright 2021 Daniel Kondor <kondor@csh.ac.at>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef CELL_IT_HELPER_HPP
#define CELL_IT_HELPER_HPP
#include <stdint.h>

template<class base>
struct it_base {
	protected:
		const base& g;
		size_t i; /* current index */
		void inc() { if(i < g.size()) i++; }
		void dec() { if(i) i--; }
		size_t add(ssize_t j) {
			if(j < 0) return sub(-1L * j);
			size_t j2 = j;
			size_t s = g.size();
			if(j2 > s || s - j2 < i) return s;
			return i + j2;
		}
		size_t sub(ssize_t j) {
			if(j < 0) return add(-1L * j);
			size_t j2 = j;
			if(j2 > i) return 0;
			return i - j2;
		}
		
		explicit it_base(const base& g_, size_t i_) : g(g_), i(i_) { }
		friend base;
		
	public:
		it_base(const it_base& it) = default;
		
		typename base::point operator *() { return g.ptid(i); }
		typename base::point operator [](ssize_t j) { return g.ptid(add(j)); }
		
		it_base& operator ++() { inc(); return *this; }
		it_base operator ++(int) { it_base tmp = *this; inc(); return tmp; }
		it_base& operator --() { dec(); return *this; }
		it_base operator --(int) { it_base tmp = *this; dec(); return tmp; }
		
		it_base& operator +=(ssize_t j) { i = add(j); return *this; }
		it_base& operator -=(ssize_t j) { i = sub(j); return *this; }
		it_base operator +(ssize_t j) { it_base tmp = *this; tmp += j; return tmp; }
		it_base operator -(ssize_t j) { it_base tmp = *this; tmp -= j; return tmp; }
		ssize_t operator -(const it_base& it) {
			if(i > it.i) return i - it.i;
			ssize_t tmp = it.i - i;
			return -1L * tmp;
		}
		
		bool operator == (const it_base& it) const { return i == it.i; }
		bool operator != (const it_base& it) const { return i != it.i; }
		bool operator < (const it_base& it) const { return i < it.i; }
		bool operator > (const it_base& it) const { return i > it.i; }
		bool operator <= (const it_base& it) const { return i <= it.i; }
		bool operator >= (const it_base& it) const { return i >= it.i; }
};

#endif

