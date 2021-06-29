/* ced.js - coherent event display
 *        - javascript functions
 *
 * $Id: ced.js 3291 2007-09-21 18:08:37Z ram $
 */

/* preload images */
function preloadImages()
{
  if (document.images)
  {
    /* preload image object */ 
    image = new Image();

    /* set images */
    url = new Array();
    url[0] = "L1_shaded_1.png";
    url[1] = "L1_shaded_05.png";
    url[2] = "H1_shaded_1.png";
    url[3] = "H1_shaded_05.png";
    url[4] = "H2_shaded_1.png";
    url[5] = "H2_shaded_05.png";
    url[6] = "V1_shaded_1.png";
    url[7] = "V1_shaded_05.png";
    url[8] = "L1_wf_noise.png";
    url[9] = "H1_wf_noise.png";
    url[10] = "H2_wf_noise.png";
    url[11] = "V1_wf_noise.png";
    url[12] = "l_tfmap_shaded_1.png";
    url[13] = "l_tfmap_shaded_05.png";
    url[14] = "l_tfmap_cluster_1.png";
    url[15] = "l_tfmap_cluster_05.png";
    url[16] = "L1_wf_strain_zoom.png";
    url[17] = "L1_wf_strain_fft.png";
    url[18] = "H1_wf_strain_zoom.png";
    url[19] = "H1_wf_strain_fft.png";
    url[20] = "H2_wf_strain_zoom.png";
    url[21] = "H2_wf_strain_fft.png";
    url[22] = "V1_wf_strain_zoom.png";
    url[23] = "V1_wf_strain_fft.png";

    /* preload images */
    var i = 0;
    for(i = 0; i < url.length; i++)
    {
      image.src = url[i];
    }
  }
}

/* function to show image */
function showImage(id, baseName, type, subType)
{
  /* update image paths */
  document.getElementById(id + "_a").href = baseName + "_" + type + "_" + subType + ".png";
  document.getElementById(id + "_img").src = baseName + "_" + type + "_" + subType + ".png";
}

/* function to show image */
function showImage1(id, baseName)
{
  /* update image paths */
  document.getElementById(id + "_a").href = baseName  + ".png";
  document.getElementById(id + "_img").src = baseName + ".png";
}

/* function to show image */
function showImage2(id, baseName, subType)
{
  /* update image paths */
  document.getElementById(id + "_a").href = baseName + "_" + subType + ".png";
  document.getElementById(id + "_img").src = baseName + "_" + subType + ".png";
}
