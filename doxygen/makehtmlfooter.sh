#!/bin/sh

# Generates the HTML footer

echo '<html>'
echo '<body>'
echo '<div id="footer" style="background-color:#E5EBF3;">'
echo '<small>'
echo '<img style="height:32px" class="footer" src="cwb_logo_circle_small_modern_alpha.png" alt="root"/></a>'
# Doxygen unconditionally adds a space in front of $DOXYGEN_ROOT_VERSION
# echo 'coherent WaveBurst'$DOXYGEN_ROOT_VERSION' - Reference Guide Generated on $datetime using 
# NOTE: $datetime must be removed to avoid regeneration of new html due only to the new date 
echo 'coherent WaveBurst'$DOXYGEN_ROOT_VERSION' - Reference Guide Generated using 
<a href="http://www.doxygen.org/index.html">
<img style="height:24px" class="footer" src="doxygen.png" alt="doxygen"/>
</a>
 '`doxygen --version`'.'
echo '</small>'
echo '</div>'
echo '</body>'
echo '</html>'
