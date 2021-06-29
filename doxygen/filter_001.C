`Begin_Macro(source)`. This parameter allows to show the macro's code in addition.
`Begin_Macro` also accept the image file type as option. "png" or "svg".
"png" is the default value. For example: `Begin_Macro(source, svg)` will show
the code of the macro and the image will be is svg format. The "width" keyword
can be added to define the width of the picture in pixel: "width=400" will
scale a picture to 400 pixel width. This allow to define large picture which
can then be scale done to have a better definition.
## In the ROOT tutorials
ROOT tutorials are also included in the ROOT documentation. The tutorials'
macros headers should look like:
~~~ {.cpp}
\file
\ingroup tutorial_hist
\notebook
Getting Contours From TH2D.
#### Image produced by `.x ContourList.C`
The contours values are drawn next to each contour.
\macro_image
#### Output produced by `.x ContourList.C`
It shows that 6 contours and 12 graphs were found.
\macro_output
#### `ContourList.C`
\macro_code
\authors  Josh de Bever, Olivier Couet
~~~
This example shows that four new directives have been implemented:
 1. `\macro_image`
 The images produced by this macro are shown. A caption can be added to document
 the pictures: `\macro_image This is a picture`
 2. `\macro_code`
 The macro code is shown.  A caption can be added: `\macro_code This is code`
 3. `\macro_output`
 The output produced by this macro is shown. A caption can be added:
 `\macro_output This the macro output`
 4. `\notebook`
   To generate the corresponding jupyter notebook. In case the tutorial does
   not generate any graphics output, the option `-nodraw` should be added.
Note that the doxygen directive `\authors` or `\author` must be the last one
of the macro header.
