=================
StereoVision Calculator
=================

.. image:: Calculator.png
   :align: center

=================
Table of contents
=================

- `Introduction`_

- `Getting started`_

- `Usage`_

  #. `Inputs`_

  #. `Results`_

  #. `Functions`_

============
Introduction
============

This program is used to calculate different parameters for a stereovision setup.
The main purpose of this program is to calculate the baseline based on a maximum depth error at a certain distance.
This however means that your software needs to have a set (estimated) disparity error. A fluctuation in its disparity error will result in a different depth error.
The reason of calculating the baseline this way is to be able to achieve an exact accuracy without the need of optimizing current algorithms.

============
Getting started
============

Clone the repository with:

.. code:: shell

    $ git clone https://github.com/ar13pit/stereovision_calculator.git

install the required python packages with:

.. code:: shell

    $ pip install -r requirements.txt

Run the bot with:


.. code:: shell

    $ python stereo_caluclator.py

That's it!

============
Usage
============

-------------------
Inputs
-------------------

| **Sensor size**
| The diagonal size of the sensor in mm or inch. (Keep in mind that sensor formats are not equal to the diagonal size, see `here <https://en.wikipedia.org/wiki/Image_sensor_format#Table_of_sensor_formats_and_sizes>`_)

| **Resolution width**
| The amount of horizontal pixels of the image

| **Resolution height**
| The amount of vertical pixels of the image

| **Focal FoV**
| The Field of View of the lens

| **Performance depth**
| The distance where you want to define a desired error

| **Performance depth error**
| The desired error that you want at the distance you previously set

| **Performance disparity**
| Select the disparity in pixels for the selected range

| **Max disparity**
| Select the maximum disparity that your software calculates (standard 128)

| **Calibration disparity error**
| Select the (expected) disparity error generated by the stereovision algorithm


-------------------
Results
-------------------

| **Focal length**
| Calculated focal length

| **Baseline**
| The distance between the two imaging sensors counted from the center of each sensor

| **Max depth**
| The maximum measurable depth (this is equal to the performance depth if the performance disparity is set to 1 pixel)

| **Depth resolution**
| The resolution of the depth image

| **Depth FoV**
| The Field of View of the depth image

-------------------
Functions
-------------------

| **Capture**
| Capture a screenshot of the program

| **Auto calculate**
| Let the program calculate the results immediately when it detects a change in the input fields

| **Calculate**
| Calculate the results

| **Plot**
| Plot the depth error chart
