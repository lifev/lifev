#! /bin/sh
#
#   SUMMARY: generate test
#     USAGE: gentest
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
# DESCRIPTION:
# ============
# Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# DESCRIP-END.
#

sed "s/test_essentialbc/$1/g" test_essentialbc_comp.qmt > $1_comp.qmt
sed "s/test_essentialbc/$1/g" test_essentialbc_exec.qmt > $1_exec.qmt
sed "s/test_essentialbc/$1/g" test_essentialbc_resu.qmt > $1_resu.qmt