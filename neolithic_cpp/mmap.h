/*  -*- C++ -*-
 * mmap.h -- simple wrapper around memory-mapped files to be used in a platform-independent way
 * 
 * Copyright 2019 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
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
 */


/* memory mapped files to read travel time data */
#ifndef MMAP_H
#define MMAP_H

#include <stdlib.h>
#include <stdint.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

class FileMapping {
	public:
		
		enum class Mode { ReadOnly, ReadWrite };
		
	
	protected:
		void* ptr; /* mapping pointer */
		size_t size; /* size of open mapping */
		off_t fsize; /* size of open file */
		int f; /* open file handle */
		Mode open_mode;
		Mode map_mode;
		bool map_shared;
		bool is_anon; /* Set if this is an anonymous mapping (i.e. memory allocation).
						 This can be useful to present a uniform interface for memory
						 that either represents a file or memory */
		FileMapping(const FileMapping&) = delete; /* this class should not be copied */
		
	public:
		FileMapping():ptr(MAP_FAILED),size(0UL),fsize(0UL),f(-1),
			open_mode(Mode::ReadOnly),map_mode(Mode::ReadOnly),map_shared(true),is_anon(false) { }
		~FileMapping() {
			close_file();
		}
		FileMapping(FileMapping&& mp):ptr(mp.ptr),size(mp.size),fsize(mp.fsize),f(mp.f),
				open_mode(mp.open_mode),map_mode(mp.map_mode),map_shared(mp.map_shared),is_anon(mp.is_anon) {
			mp.ptr = MAP_FAILED;
			mp.f = -1;
			mp.size = 0UL;
			mp.fsize = 0UL;
		}
		
		/**
		 * Open an existing file, name given by fn.
		 * The parameter mode gives the file access mode.
		 * If create is true, a new file is created if it does not already exists;
		 * otherwise it is an error if the file does not exist.
		 * Note: the combination mode == Mode::ReadOnly and create == true makes
		 * little sense, but is allowed.
		 * 
		 * Any existing mapping is discarded and any open file is closed. */
		bool open_file(const char* fn, Mode mode = Mode::ReadOnly, bool create = false) {
			close_file();
			open_mode = mode;
			mode_t create_mode = 00644;
			int flags = O_CLOEXEC | (mode == Mode::ReadWrite ? O_RDWR : O_RDONLY);
			if(create) flags |= O_CREAT;
			f = open(fn,flags,create_mode);
			if(f == -1) return false;
			
			/* check and store file size */
			struct stat sb;
			if(fstat(f,&sb)) {
				close_file();
				return false;
			}
			fsize = sb.st_size;
			return true;
		}
		
		/**
		 * Sets this instance to use an anonymous mapping (i.e. allocating memory from the OS).
		 * Use change_file_size() to set the size of memory requested before calling map_file().
		 * Note: contents of a mapping created this way are lost when it is unmapped or closed.
		 * Any existing mapping is discarded and any open file is closed. */
		bool open_anon() {
			close_file();
			is_anon = true;
			return true;
		}
		
		bool change_file_size(off_t new_size) {
			/* cannot change size if no file is open */
			if(f == -1) {
				if(is_anon) { fsize = new_size; return true; }
				else return false;
			}
			/* cannot change size if we have no write access */
			if(open_mode != Mode::ReadWrite) return false;
			/* should not change size if we have an open mapping */
			if(ptr != MAP_FAILED) return false;
			if(ftruncate(f, new_size)) return false;
			fsize = new_size;
			return true;
		}
		
		bool map_file(Mode mode = Mode::ReadOnly, bool shared = true) {
			if(f == -1 && !is_anon) return false;
			if(mode != map_mode || shared != map_shared) unmap_file();
			map_mode = mode;
			map_shared = (shared && !is_anon);
			
			int prot = PROT_READ;
			if(mode == Mode::ReadWrite) prot |= PROT_WRITE;
			
			int flags = shared ? MAP_SHARED : MAP_PRIVATE;
			if(is_anon) flags |= MAP_ANONYMOUS;
			
			if(ptr == MAP_FAILED) ptr = mmap(0,fsize,prot,flags,f,0);
			if(ptr != MAP_FAILED) {
				size = fsize;
				return true;
			}
			else {
				size = 0UL;
				return false;
			}
		}
		
		/* Try to ensure that all data is read into memory. For anonymous mappings, this will likely result in
		 * physical memory being allocated. Note: this function assumes a page size of 4 KiB.
		 * Also, we assume that sizeof(char) == 1. */
		void prefault() {
			if(!is_mapped()) return;
			const size_t page_size = 4096;
			volatile uint64_t dummy = 0;
			size_t max = size;
			if(max > std::numeric_limits<size_t>::max() - page_size) max -= page_size;
			const char* base = (const char*)ptr;
			for(size_t i = 0; i < max; i += page_size) {
				dummy += base[i];
			}
		}
		
		void unmap_file() {
			if(ptr != MAP_FAILED) munmap(ptr,size);
			ptr = MAP_FAILED;
			size = 0UL;
		}
		
		void close_file() {
			unmap_file();
			if(f != -1) close(f);
			f = -1;
		}
		
		bool set_mem_mode(Mode mode) {
			if(mode == map_mode) return true;
			if(ptr == MAP_FAILED) return false;
			int prot = PROT_READ;
			if(mode == Mode::ReadWrite) prot |= PROT_WRITE;
			if(mprotect(ptr, size, prot)) return false;
			map_mode = mode;
			return true;
		}
		
		off_t file_size() const { return fsize; }
		const void* get_mapping() const { return (ptr == MAP_FAILED) ? nullptr : ptr; }
		void* get_mapping() { return (ptr == MAP_FAILED) ? nullptr : ptr; }
		bool is_mapped() const { return (ptr != MAP_FAILED); }
		bool is_open() const { return (f != -1); }
		bool is_anon_map() const { return is_anon; }
		Mode file_mode() const { return open_mode; }
		Mode mem_mode() const { return map_mode; }
		bool is_shared() const { return map_shared; }
};


/* templated mapping for specific data type */
template<class T>
class FileMappingT : public FileMapping {
	public:
		size_t size() const { return fsize / sizeof(T); }
		T& operator [] (size_t i) {
			T* tmp = (T*)get_mapping();
			return tmp[i];
		}
		const T& operator [] (size_t i) const {
			T* tmp = (T*)get_mapping();
			return tmp[i];
		}
		T* data() { return (T*)get_mapping(); }
		const T* data() const { return (const T*)get_mapping(); }
};


#endif

