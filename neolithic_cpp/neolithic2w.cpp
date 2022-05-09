/*
 * neolithic2w.cpp -- first version of neolithic simulation, abstract space with warfare
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

#include "neolithic_common.hpp"

int main(int argc, char **argv)
{
	neolithic_sim<true> sim;
	sim.run(argc, argv);
}

