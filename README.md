# MATH5610_TermProject

When running/piping the documents in Unix it is CRUCIAL that you run/pipe the <satellite.py> and <receiver.py> with the prefix "python3" specified as the language and NOT just "python". Running the classes with the "python" prefix will yield inaccurate results. We suspect there was an update in python3 that modified a built-in method, hence the lack of backwards compatibility.

    cat bt0.dat | java vehicle | python3 satellite.py | python3 receiver.py

Note: The NumPy library was used so it must be installed on the machine running the program. The class files will automatically import the NumPy library but it MUST be installed for the scripts to run.

We tested and developed our code on David's personal computer and periodically ran it through the school's Unix system to ensure compatibility.

                    c-jsdk@xserver.math.utah.edu

To help ensure accurate results we used <satellite_test.py> and <receiver_test.py>. The two programs simply just compare two files and compare the locations and times. If locations on both files are within a centimeter, the test passes. If the times or locations differ too much, the test fails. You can test the two files by typing the file names with a space inbetween, like below.

                   <this_file>.log <other_file>.log




