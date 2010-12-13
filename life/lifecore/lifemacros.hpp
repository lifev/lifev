//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
/*!
    @file
    @brief Base utilities operating on meshes
    @date 2005-01-24

    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @contributor
    @maintainer Simone Deparis <simone.deparis@epfl.ch>

    @date 08-02-2004

    @subsection hints LifeV C++ Compiler Hints

    -# #INLINE
    -# #LIFEV_RESTRICT

    @subsection attribute_macro LifeV Attribute Macros

    -# #LIFEV_EXPORT and #LIFEV_NO_EXPORT
    -# #LIFEV_PACKED
    -# #LIFEV_DEPRECATED
    -# #LIFEV_ISLIKELY and #LIFEV_ISUNLIKELY

 */

#ifndef LIFEMACROS_HPP
#define LIFEMACROS_HPP 1

/**
   @def LIFEV_CONSTRUCTOR_BEGIN(Area,x)
   Inform that the constructor of the class x has started
 */
#define LIFEV_CONSTRUCTOR_BEGIN(Area, A) Debug( Area ) << "Constructor of " << A << " begins\n";
#define LIFEV_CONSTRUCTOR(Area,A) LIFEV_CONSTRUCTOR_BEGIN(Area,A)
#define CONSTRUCTOR(A) LIFEV_CONSTRUCTOR_BEGIN(20000,A)

/**
   @def LIFEV_CONSTRUCTOR_END(Area,x)
   Inform that the constructor of the class x has ended
 */
#define LIFEV_CONSTRUCTOR_END(Area,A) Debug( Area ) << "Constructor of " << A << " ends\n";

/**
   @def LIFEV_DESTRUCTOR_BEGIN(Area,x)
   Inform that the destructor of the class x has started
 */
#define LIFEV_DESTRUCTOR_BEGIN(Area,A) Debug( Area ) << "Destructor of " << A << " begins\n";
#define LIFEV_DESTRUCTOR(Area,A) LIFEV_DESTRUCTOR_END(Area,A)
#define DESTRUCTOR(A) LIFEV_DESTRUCTOR_BEGIN(20000,A)

/**
   @def LIFEV_DESTRUCTOR_END(Area,x)
   Inform that the destructor of the class x has started
 */
#define LIFEV_DESTRUCTOR_END(Area,A) Debug( Area ) << "Destructor of " << A << " ends\n";


/**
   @def INLINE

   Alias to the C/C++ keyword \c inline
 */
#define INLINE inline

/**
   @def LIFEV_RESTRICT
   @brief C99 feature of restricted(not aliased) pointers and references

   As with gcc, g++ understands the C99 feature of restricted
   pointers, specified with the \c __restrict__, or __restrict type
   qualifier. Because you cannot compile C++ by specifying the
   -std=c99 language flag, restrict is not a keyword in C++.


   In addition to allowing restricted pointers, you can specify
   restricted references, which indicate that the reference is not
   aliased in the local context.

   @code
   void fn (int *__restrict__ rptr, int &__restrict__ rref)
   {
   ...
   }
   @endcode

   In the body of \c fn, \c rptr points to an unaliased integer and \c rref
   refers to a (different) unaliased integer.


   You may also specify whether a member function's this pointer is
   unaliased by using \c __restrict__ as a member function qualifier.

   @code
   void T::fn () __restrict__
   {
   ...
   }
   @endcode

   Within the body of T::fn, this will have the effective definition
   <tt>T* __restrict__</tt> const this. Notice that the interpretation of a
   \c __restrict__ member function qualifier is different to that of
   const or volatile qualifier, in that it is applied to the pointer
   rather than the object. This is consistent with other compilers
   which implement restricted pointers.

   As with all outermost parameter qualifiers, \c __restrict__ is ignored
   in function definition matching. This means you only need to
   specify \c __restrict__ in a function definition, rather than in a
   function prototype as well.

   In order to ensure that the code is portable to other compiler than
   gcc/g++ a macro has been defined LIFEV_RESTRICT that is equal to
   __restrict__ if the compiler supports it.
 */
#define LIFEV_RESTRICT __restrict__





/**
   @def LIFEV_EXPORT
   @brief Load time improvements for DSO libraries

   Here are a few explanations why this is useful.  For more info
   checkout http://www.nedprod.com/programs/gccvisibility.html

   -# It very substantially improves load times of your DSO (Dynamic
   Shared Object) For example, the TnFOX Boost.Python bindings library
   now loads in eight seconds rather than over six minutes!

   -# It lets the optimiser produce better code PLT indirections (when a
   function call or variable access must be looked up via the Global
   Offset Table such as in PIC code) can be completely avoided, thus
   substantially avoiding pipeline stalls on modern processors and thus
   much faster code. Furthermore when most of the symbols are bound
   locally, they can be safely elided (removed) completely through the
   entire DSO. This gives greater latitude especially to the inliner
   which no longer needs to keep an entry point around "just in case".

   -# It reduces the size of your DSO by 5-20% ELF's exported symbol
   table format is quite a space hog, giving the complete mangled symbol
   name which with heavy template usage can average around 1000
   bytes. C++ templates spew out a huge amount of symbols and a typical
   C++ library can easily surpass 30,000 symbols which is around 5-6Mb!
   Therefore if you cut out the 60-80% of unnecessary symbols, your DSO
   can be megabytes smaller!

   -# Much lower chance of symbol collision The old woe of two libraries
   internally using the same symbol for different things is finally
   behind us with this patch. Hallelujah!

   here is an example on how to use them
   @code
   int LIFEV_NO_EXPORT foo;
   int LIFEV_EXPORT bar;

   extern "C" LIFEV_EXPORT void function(int a);

   class LIFEV_EXPORT SomeClass
   {
     int c;

     // Only for use within this DSO
     LIFEV_NO_EXPORT void privateMethod();

    public:

     Person(int _c) : c(_c) { }
     static void foo(int a);
    };
   @endcode
*/
/**
   @def LIFEV_NO_EXPORT

   Counterpart to #LIFEV_EXPORT.
 */
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 > 2)
#define LIFEV_EXPORT __attribute__ ((visibility("default")))

#define LIFEV_NO_EXPORT __attribute__ ((visibility("hidden")))
#else
#define LIFEV_EXPORT
#define LIFEV_NO_EXPORT
#endif

/**
   @def LIFEV_PACKED
   The LIFEV_PACKED can be used to hint the compiler that a particular
   structure or class should not contain unnecessary paddings.

   Here is an explanation from http://sig9.com/articles/gcc-packed-structures

   GCC allows you to specify attributes of variables and structures
   using the keyword \c __attribute__, the syntax of which is
   \c __attribute__((attribute list)). One such attribute is \c __packed__
   which specifies that

   a variable or structure field should have the smallest possible
   alignment--one byte for a variable, and one bit for a field, unless
   you specify a larger value with the aligned attribute.


   which means that GCC will not add any of the zero's for padding (for
   memory alignement) and make variables or fields immediately next to
   each other. For example, here are some things I tried out -- I created
   a C source file - \c test.c

   @code
   struct test_t {
   int  a;
   char b;
   int  c;
   } ;

   struct test_t test = { 10, 20, 30};
   @endcode

   And compiled it with the -S option (ie to generate the assembly
   equivalent of the code generated).

   @code
      .file "t.cpp"
      .globl test
        .data
        .align 4
        .type test, @object
        .size test, 12
      <b>test:
      .long 10
      .byte 20
      .zero 3
      .long 30</b>
      .section .note.GNU-stack,"",@progbits
      .ident   "GCC: (GNU) 3.3.5 (Debian 1:3.3.5-6)"
   @endcode

   Notice the emphasized code. You can see that the structure "test"
   is being declared. First the field "a" (int) as .long 10 followed
   by "b" (char) as .byte 20. To keep the fields' word alignment,
   notice that GCC has added 3 zero bytes (.zero 3) before field "c"
   (int) which is declared as .long 30. This makes the effective
   sizeof struct test_t as 12 instead of the expected 9. Then I tried
   with the __packed__ attribute -

   @code
   struct test_t {
   int  a;
   char b;
   int  c;
   } LIFEV_PACKED

   struct test_t test = { 10, 20, 30};
   @endcode

   and the "-S" output I got after compiling was

   @code
   .file "t.cpp"
   .globl test
     .data
     .type test, @object
     .size test, 9
   test:
     .long 10
     .byte 20
     .long 30
   .section .note.GNU-stack,"",@progbits
   .ident   "GCC: (GNU) 3.3.5 (Debian 1:3.3.5-6)"
   @endcode

   in which the zeros are missing making the sizeof structure test_t =
   9. Always remember that memory alignment is *good* even if it
   compromises space, so think twice before using this attribute. It
   is generally useful when you want to assign a structure to a block
   of memory and manipulate it through the fields of a structure.
 */
#ifdef __GNUC__
#define LIFEV_PACKED __attribute__((__packed__))
#else
#define LIFEV_PACKED
#endif

/**
   The LIFEV_DEPRECATED macro can be used to trigger compile-time warnings
   with gcc >= 3.2 when deprecated functions are used.

   For non-inline functions, the macro gets inserted at the very end of the
   function declaration, right before the semicolon:

   @code
   DeprecatedConstructor() LIFEV_DEPRECATED;
   void deprecatedFunctionA() LIFEV_DEPRECATED;
   int deprecatedFunctionB() const LIFEV_DEPRECATED;
   @endcode

   Functions which are implemented inline are handled differently: for them,
   the LIFEV_DEPRECATED macro is inserted at the front, right before the return
   type, but after "static" or "virtual":

   @code
   LIFEV_DEPRECATED void deprecatedInlineFunctionA() { .. }
   virtual LIFEV_DEPRECATED int deprecatedInlineFunctionB() { .. }
   static LIFEV_DEPRECATED bool deprecatedInlineFunctionC() { .. }
   @end

   You can also mark whole structs or classes as deprecated, by inserting the
   LIFEV_DEPRECATED macro after the struct/class keyword, but before the
   name of the struct/class:

   @code
   class LIFEV_DEPRECATED DeprecatedClass { };
   struct LIFEV_DEPRECATED DeprecatedStruct { };
   @endcode
*/
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 >= 2)
# define LIFEV_DEPRECATED __attribute__ ((deprecated))
#else
# define LIFEV_DEPRECATED
#endif

/**
   @def LIFEV_ISLIKELY(x)
   The LIFEV_ISLIKELY macro tags a boolean expression as likely to evaluate to
   'true'. When used in an if ( ) statement, it gives a hint to the compiler
   that the following codeblock is likely to get executed. Providing this
   information helps the compiler to optimize the code for better performance.
   Using the macro has an insignificant code size or runtime memory footprint impact.
   The code semantics is not affected.

   @note
   Providing wrong information ( like marking a condition that almost never
   passes as 'likely' ) will cause a significant runtime slowdown. Therefore only
   use it for cases where you can be sure about the odds of the expression to pass
   in all cases ( independent from e.g. user configuration ).

   @par
   The LIFEV_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.

   @note
   Do NOT use ( !LIFEV_ISLIKELY(foo) ) as an replacement for LIFEV_ISUNLIKELY !

   @code
   if ( LIFEV_ISUNLIKELY( testsomething() ) )
       abort();     // assume its unlikely that the application aborts
   @endcode
*/
/**
   @def LIFEV_ISUNLIKELY(x)
   Counterpart to #LIFEV_ISLIKELY
   The LIFEV_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.
 */
#if __GNUC__ - 0 >= 3
# define LIFEV_ISLIKELY( x )    __builtin_expect(!!(x),1)
# define LIFEV_ISUNLIKELY( x )  __builtin_expect(!!(x),0)
#else
# define LIFEV_ISLIKELY( x )   ( x )
# define LIFEV_ISUNLIKELY( x )  ( x )
#endif

#endif /* LIFEMACROS_HPP */
