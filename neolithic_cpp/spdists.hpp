/*
 * spdists.hpp -- store a (symmetric) matrix of distances among a set of entities
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


#ifndef SPDISTS_HPP
#define SPDISTS_HPP
#include <memory>
#include <stdint.h>

template<class dtype>
class spdists {
	protected:
		const dtype* d = nullptr; /* distance matrix (read-only) */
		std::shared_ptr<dtype[]> pd; /* read-write pointer to the same as above (only set if distances were allocated by this class) */
		size_t size1 = 0; /* size of the matrix (number of elements); size of d is size*(size+1) / 2 */
		
		/* get the index where the distance between node i and j is stored in the matrix
		 * (throws exception if i or j is above size); */
		size_t get_ix(size_t i, size_t j) const {
			if(i >= size1 || j >= size1) throw std::runtime_error("spdists::get_ix(): node indeces out of range!\n");
			if(i < j) {
				size_t tmp = i;
				i = j;
				j = tmp;
			}
			/* here i >= j */
			return (i*(i+1UL)) / 2UL + j;
		}
		
		/* set the distance between two nodes -- this assumes that pd is set
		 * (this is only called internally) */
		void set_dist0(size_t i, size_t j, const dtype& x) {
			pd[get_ix(i, j)] = x;
		}
		
		constexpr static uint64_t header_base = 0x635b054329591700UL;
		constexpr static uint64_t get_header_dtype() {
			throw std::runtime_error("Unsupported distance data type!");
			return (uint64_t)-1; /* hack for clang */
		}
	
	public:
		size_t size() const { return size1; }
		bool has_data() const { return (d != nullptr) && size1 > 0UL; }
		const dtype& get_dist(size_t i, size_t j) const {
			return d[get_ix(i, j)];
		}
		
		/* allocate a new distance matrix for a given size, with the
		 * default element set to the second argument */
		void allocate(size_t s, const dtype& default_dist = dtype()) {
			size1 = s;
			size_t tmp = (size1 * (size1 + 1UL)) / 2UL;
			pd.reset(new dtype[tmp]);
			d = pd.get();
			for(size_t i = 0; i < tmp; i++) pd[i] = default_dist;
		}
		
		/* free up the stored distribution */
		void reset() {
			pd.reset();
			d = nullptr;
			size1 = 0;
		}
		
		/* set the distance between two nodes */
		void set_dist(size_t i, size_t j, const dtype& x) {
			if(!pd) throw std::runtime_error("spdists::set_dist(): storage needs to be allocated first!\n");
			set_dist0(i, j, x);
		}
		
		/* write this matrix to a binary file
		 * returns the number of bytes written or 0 on error */
		size_t write(FILE* f) const {
			if(!(size1 && d)) throw std::runtime_error("spdists::write(): nothing to write!\n");
			size_t ret = 0;
			uint64_t tmp1 = header_base + get_header_dtype();
			if(fwrite(&tmp1, sizeof(uint64_t), 1, f) != 1) return 0;
			ret += sizeof(uint64_t);
			if(fwrite(&size1, sizeof(size_t), 1, f) != 1) return 0;
			ret += sizeof(size_t);
			
			/* potential padding */
			size_t align = alignof(dtype);
			if(ret % align) {
				uint8_t padding[align];
				for(size_t i = 0; i < align; i++) padding[i] = 0;
				size_t pad = align - (ret % align);
				if(fwrite(padding, 1, pad, f) != pad) return 0;
				ret += pad;
			}
			
			tmp1 = (size1 * (size1 + 1UL)) / 2UL;
			if(fwrite(d, sizeof(dtype), tmp1, f) != tmp1) return 0;
			ret += tmp1 * sizeof(dtype);
			return ret;
		}
		
		/* set up the distance matrix from the given memory area
		 * (corresponding to a file written by write() and opened by mmap());
		 * returns the number of bytes consumed
		 * note: it is the caller's responsibility to ensure that the
		 * given memory stays valid while this class is used */
		size_t read(const void* base, size_t base_size) {
			d = nullptr;
			size1 = 0;
			pd.reset();
			size_t ret = 0;
			
			if(base_size < 2*sizeof(size_t)) throw std::runtime_error("spdists::read(): input is too small!\n");
			uint64_t header = *((const uint64_t*)base);
			if(header != header_base + get_header_dtype()) throw std::runtime_error("spdists::read(): invalid header!\n");
			const uint8_t* base2 = (const uint8_t*)base;
			ret += sizeof(uint64_t);
			size1 = *((const size_t*)(base2 + ret));
			ret += sizeof(size_t);
			
			/* potential padding */
			size_t align = alignof(dtype);
			if(ret % align) {
				size_t pad = align - (ret % align);
				ret += pad;
			}
			
			size_t tmp = (size1 * (size1 + 1UL)) / 2UL;
			tmp *= sizeof(dtype);
			if(tmp + ret > base_size) {
				size1 = 0;
				throw std::runtime_error("spdists::read(): input is too small!\n");
			}
			d = (const dtype*)(base2 + ret);
			return ret + tmp;
		}
};


/* header bits identifying possible data types */
template<> constexpr uint64_t spdists<int16_t>::get_header_dtype()  { return 1UL; }
template<> constexpr uint64_t spdists<uint16_t>::get_header_dtype() { return 2UL; }
template<> constexpr uint64_t spdists<int32_t>::get_header_dtype()  { return 3UL; }
template<> constexpr uint64_t spdists<uint32_t>::get_header_dtype() { return 4UL; }
template<> constexpr uint64_t spdists<int64_t>::get_header_dtype()  { return 5UL; }
template<> constexpr uint64_t spdists<uint64_t>::get_header_dtype() { return 6UL; }
template<> constexpr uint64_t spdists<float>::get_header_dtype()    { return 7UL; }
template<> constexpr uint64_t spdists<double>::get_header_dtype()   { return 8UL; }


#endif

