#!/bin/bash
#
# -----------------------------------------------------------------------
# This file is part of the hydrangea tools package
# Copyright (C) 2019 Yannick Bahe (bahe@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------

gcc -fPIC -g -O3 -c -Wall main.c 
gcc -shared -lm -lc -o sumbins.so main.o
cp sumbins.so ../../sumbins.so

