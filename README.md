<h1><i>Mycoplasma genitalium</i> whole-cell network</h1>

In order to make the whole-cell network of Mycoplasma genitalium reproducible, we made available all the scripts developed to build the network. There are two Python 2.7 scripts, one to read the database and create the model in SBML format and other to convert the SBML file to a GML file. The usage of each one is described below.

<b>build_mg_wcnetwork.py</b> This script require previous installation of Python 2.7 compiler, SBML API lib and MySQLdb lib for Python 2.7 and the EMBOSS transeq software. A MySQL database also must be running with the \textit{Mycoplasma genitalium} Knowledge Database. This database can be obtained at https://simtk.org/frs/download.php?file_id=3426 with the file \verb|data.sql| inside the zip file. To run the script, use the command "python2.7 build_mg_wcnetwork.py.py -h" and the instructions to run will be shown.

<b>sbml2gml.py</b> This script require previous installation of Python 2.7 compiler, SBML API lib and networkx lib for Python 2.7. To use the script, use the command "python2.7 sbml2gml.py -h" and the instructions to run will be shown.
