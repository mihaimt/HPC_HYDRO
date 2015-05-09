Instructions for visualizing the HYDRO output
============================================
=========================================== 


Plots:
============================================
python output_to_plots.py output_name

[based on Rafael parse_vtk.py] 

Files generated:
[output_name]_d.png
[output_name]_p.png
[output_name]_u.png
[output_name]_v.png
===========================================


Movies:
==============================================
python output_to_movies.py output_list.txt 100
				|           |
				|	   time delay between the movie frames
				|
For movies you will need an output_list.txt file with the structure:
----------------------------------------------------------
out1
out2
out3
....
outn 
-----------------------------------------------------------
Example: output_list.txt
------------------------------------------------
outputvtk_00006.vts
outputvtk_00001.vts
outputvtk_00007.vts
outputvtk_00002.vts
outputvtk_00008.vts
outputvtk_00003.vts
outputvtk_00009.vts
outputvtk_00004.vts
outputvtk_00010.vts
outputvtk_00005.vts
===============================================


Order does not matter since the movie creator reorders according to the number
found in the name of the file.

Files generated:
[output_name]_d.mpg
[output_name]_p.mpg
[output_name]_u.mpg
[output_name]_v.mpg
===========================================

Sugestion:
Use command

xdg-open movie_name.mpg

for proper viewing of the movie.

