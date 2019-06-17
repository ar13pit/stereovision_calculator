from __future__ import (print_function, division, absolute_import,
                        unicode_literals)
from builtins import *
import tkinter as tk
from tkinter import ttk
from ttkthemes import ThemedTk
from pyscreenshot import grab
import matplotlib.pyplot as plt
from matplotlib import style
import math


def baselineCalculator(focal_length, disparity, depth):
    """
    Function to calculate the baseline (in m) of the stereo camera given:
    :param: disparity: Difference between the x pixels of both images (in px)
    :param: depth: The depth corresponding to the disparity (in m)
    :param: focal_length: The focal length of the system (in px)
    """
    return (disparity * depth / focal_length)

def disparityToDepth(baseline, focal_length, disparity):
    """
    Function to calculate the Depth (in m) of a image point given:
    :param: baseline: The distance between the camera centers (in m)
    :param: focal_length: The focal length of the system (in px)
    :param: disparity: Difference between the x pixels of both images (in px)
    """
    return (baseline * focal_length /disparity)

def depthToDisparity(baseline, focal_length, depth):
    """
    Function to calculate the disparity between image points (in px) given:
    :param: baseline: The distance between the camera centers (in m)
    :param: focal_length: The focal length of the system (in px)
    :param: depth: Depth of the image point (in m)
    """
    return (baseline * focal_length /depth)

class StereoVisionCalculator(object):
    def __init__(self):
        """
        Class constructor of StereoVisionCalculator.
        This method initalizes the tkinter GUI and the object variables to
        store the user input and calculation results.
        """
        # Focal length calculator parameters
        self.sensor_size = None
        self.img_width = None
        self.img_height = None
        self.focal_fov = None
        self.focal_length = None

        # Design limits
        self.perf_baseline_min = None
        self.perf_baseline_max = None
        self.perf_depth_max = None
        self.perf_depth_min = None
        self.perf_depth_error = None
        self.perf_disp_max = None
        self.perf_disp_min = None
        self.perf_disp_calibration_error = None

        # Baseline calculator results
        self._max_depth = None
        self._min_depth = None
        #self.max_depth_error = None
        #self.max_disparity = None
        #self.calibration_disparity_error = None
        #self.min_disparity = None
        self._baseline = None
        #self.min_depth = None

        # Initialize the complete GUI
        self._initializeGUI()

    def _initializeGUI(self):
        """
        Method to setup the StereoVision Calculator GUI using tkinter
        """
        self.root = ThemedTk(theme="arc")
        self.root.tk_setPalette(background='#f5f6f7')
        self.root.title("StereoVision Calculator")
        self.root.resizable(0, 0)  # Don't allow resize

        self.check = [tk.IntVar(self.root) for _ in range(2)]  # Checkboxes
        self.pmenu = [tk.StringVar(self.root) for _ in range(7)]  # Popup menus
        self.results = [tk.DoubleVar(self.root) for _ in range(6)]  # Results
        self.entries = [ttk.Entry(self.root) for _ in range(12)]  # Entries
        self.entries[9]["state"] = "disabled"
        self.entries[9]["width"] = 14

        self.results[4].set('0 x 0')
        self.results[5].set('0° x 0°')

        # Entries
        for row, entry in enumerate(self.entries, start=1):
            entry.grid(row=row, column=1, sticky="W")

        labels = {
            'Sensor size': '',
            'Resolution width': 'pixels',
            'Resolution height': 'pixels',
            'Focal FoV': '',
            'Min baseline': '',
            'Max baseline': '',
            'Performance min depth': '',
            'Performance depth': '',
            'Performance depth error': '',
            'Performance disparity': 'pixels',
            'Max disparity': 'pixels',
            'Calibration disparity error': 'pixel',
            'Focal length': 'mm',
            'Baseline': 'mm',
            'Min depth': 'cm',
            'Max depth': 'm',
            'Depth resolution': 'pixels',
            'Depth FoV': 'deg'
        }

        for x, (val, key) in enumerate(labels.items()):
            ttk.Label(self.root, text=val).grid(
                row=x + 1, sticky="W", pady=5)
            if key:
                ttk.Label(self.root, text=key).grid(
                    row=x + 1, column=2, sticky="W")

        # Buttons
        ttk.Checkbutton(self.root,
                        text="Auto calculate",
                        variable=self.check[0],
                        command=self._callback).grid(row=0, sticky="W", pady=5)
        ttk.Checkbutton(self.root,
                        text="",
                        variable=self.check[1],
                        command=self._disp).grid(row=10, column=1, sticky="E")
        ttk.Button(self.root,
                   text="Capture",
                   width=12,
                   command=self._capture).grid(row=0, column=2, sticky="W")
        ttk.Button(self.root,
                   text="Calculate",
                   width=12,
                   command=self._callback).grid(row=19, sticky="W")
        ttk.Button(self.root,
                   text="Plot",
                   width=12,
                   command=self._plot).grid(row=19, column=2, sticky="W")

        # Dropdown menus
        choices1 = {'', 'mm', 'in'}
        choices2 = {'', 'Horizontal', 'Vertical', 'Diagonal'}
        choices3 = {'', 'mm', 'cm', 'm', 'in', 'ft'}
        self.popupMenu = []
        for x in range(0, 7):
            if x == 0:
                choices = choices1
                self.pmenu[x].set('mm')
                a = 1
            elif x == 1:
                choices = choices2
                self.pmenu[x].set('Horizontal')
                a = 3
            else:
                choices = choices3
                self.pmenu[x].set('mm')
            self.popupMenu.append(ttk.OptionMenu(self.root, self.pmenu[x], *choices))
            self.popupMenu[x].grid(row=x + a, column=2, sticky="W")

        # Results
        for x in range(0, len(self.results)):
            ttk.Label(self.root, textvariable=self.results[x]).grid(
                row=x + 13, column=1, sticky="W")

        col_count, row_count = self.root.grid_size()

        for col in range(col_count):
            self.root.grid_columnconfigure(col, pad=2)

    def mainloop(self):
        self.root.mainloop()

    def calculateToMeter(self, variable, conversion_from):
        if conversion_from is "mm":
            result = variable / 1000
        elif conversion_from is "cm":
            result = variable / 100
        elif conversion_from is "m":
            result = variable
        elif conversion_from is "in":
            result = variable * 0.0254
        elif conversion_from is "ft":
            result = variable * 0.3048

        return result

    def _focalLengthCalculator(self, sensor_size, size, img_width, img_height,
                               focal_fov, fov_type):
        """
        Function to calculate the focal length of the imaging sensor given:
        :param: sensor_size: The diagonal sensor size
        :param: size: The measurement system of the sensor (metric/imperial)
        :param: img_width: The amount of pixels in width
        :param: img_height:  The amount of pixels in height
        :param: focal_fov: The field of view of the lens
        :param: fov_type: Horizontal, vertical or diagonal FoV
        """
        if size == 'inch':
            sensor_size = sensor_size * 25.4
        else:
            sensor_size = sensor_size

        ratio = img_height / img_width
        sensor_width_mm = math.sqrt(sensor_size * sensor_size /
                                    (1.0 + ratio * ratio))
        sensor_height_mm = math.sqrt(sensor_size * sensor_size /
                                     (1.0 + 1.0 / (ratio * ratio)))

        roi_width_mm = sensor_width_mm  # * roi_width / img_width
        roi_height_mm = sensor_height_mm  # * roi_height / img_height
        roi_diagonal_mm = math.sqrt(roi_height_mm * roi_height_mm +
                                    roi_width_mm * roi_width_mm)

        fov = focal_fov / 180 * math.pi
        atanInner = math.tan(fov * 0.5)

        try:
            if fov_type == 'Horizontal':
                f_mm = roi_width_mm / (2 * atanInner)
            elif fov_type == 'Vertical':
                f_mm = roi_height_mm / (2 * atanInner)
            elif fov_type == 'Diagonal':
                f_mm = roi_diagonal_mm / (2 * atanInner)

            pixel_size_mm = roi_width_mm / img_width
            f_pixel = f_mm / pixel_size_mm
        except ZeroDivisionError:
            f_mm = 0
            f_pixel = 0

        return f_mm, f_pixel, roi_width_mm, roi_height_mm

    def _baselineCalculator(self, focal_length, max_depth, max_depth_error,
                            max_disparity, calibration_disparity_error,
                            min_disparity):
        """
        Function to calculate the baseline of the stereo camera given:
        :param: min_disparity_measured: The least measurable disparity of the
            system
        :param: max_depth: The real depth corresponding to the
            min_disparity_measured
        :param: max_depth_error: The maximum allowed error in the measurement of
            depth at max_depth
        :param: disparity_range: Total number of disparity points allowed. So
            maximum disparity becomes min_disparity_measured + disparity_range
        :param: focal_length: The focal length of the system in terms of pixels
        """
        baseline = ((min_disparity + calibration_disparity_error) *
                    (max_depth + max_depth_error)) / focal_length

        min_depth = (baseline * focal_length) / max_disparity

        return baseline, min_depth

    def _depthErrorCalculator(self, depth, baseline, focal_length,
                              max_disparity_error):
        """
        Method to calculate the max_depth_error for a given depth for
        pre-determined baseline and focal length
        :param: depth
        :return: max_depth_error
        """
        disparity_real = baseline * focal_length / depth
        disparity_measured = disparity_real + max_disparity_error

        depth_measured = baseline * focal_length / disparity_measured
        return abs(depth_measured - depth)

    def _depthCalculator(self, roi_width, roi_height, roi_width_mm, roi_height_mm, img_width, d_max, f_mm):
        roi_full_pixel = str(roi_width - d_max) + ' x ' + str(roi_height)
        h_full_angle = 2 * math.atan((1 * roi_width_mm * (img_width - d_max) / img_width) / (2 * f_mm))
        v_angle = 2 * math.atan(roi_height_mm / (2 * f_mm))

        FoV_h = round(h_full_angle / math.pi * 180, 1)
        FoV_v = round(v_angle / math.pi * 180, 1)
        FoV = str(FoV_h) + '° x ' + str(FoV_v) + '°'
        return FoV, roi_full_pixel

    def _callback(self):

        if (self.entries[0].get() and self.entries[1].get() and
                self.entries[2].get() and self.entries[3].get()):

            sensor_size = float(self.entries[0].get())
            size = self.pmenu[0].get()
            img_width = int(self.entries[1].get())
            img_height = int(self.entries[2].get())
            focal_fov = float(self.entries[3].get())
            fov_type = self.pmenu[1].get()

            f_mm, f_pixel, roi_width_mm, roi_height_mm = self._focalLengthCalculator(
                sensor_size, size, img_width, img_height, focal_fov, fov_type)
            self.results[0].set(round(f_mm, 2))

            if (self.entries[4].get() and self.entries[5].get() and
                    self.entries[6].get() and self.entries[7].get()):
                focal_length = f_pixel
                min_baseline = calculateToMeter(
                    float(self.entries[4].get()), self.pmenu[2].get)
                max_baseline = calculateToMeter(
                    float(self.entries[5].get()), self.pmenu[3].get)
                min_depth = calculateToMeter(
                    float(self.entries[6].get()), self.pmenu[6].get)
                perf_depth = calculateToMeter(
                    float(self.entries[7].get()), self.pmenu[7].get)
                perf_depth_err = calculateToMeter(
                    float(self.entries[8].get()), self.pmenu[8].get)
                perf_disp = 1 if not self.entries[9].get() else int(
                    self.entries[9].get())
                max_disp = int(self.entries[10].get())
                calibration = float(self.entries[11].get())

                # baseline, min_depth, max_depth = self._baselineCalculator(
                #     focal_length, max_depth, max_depth_error,
                #     max_disparity, calibration_disparity_error,
                #     min_disparity)

                self.results[1].set(round(baseline * 1000, 2))
                self.results[2].set(round(min_depth * 100, 2))
                self.results[3].set(round(max_depth * 1000, 2))

            if self.entries[10].get():
                d_max = int(self.entries[10].get()) - 1
                depth_fov, depth_res = self._depthCalculator(
                    img_width, img_height, roi_width_mm, roi_height_mm, img_width, d_max, f_mm)
                self.results[4].set(depth_res)
                self.results[5].set(depth_fov)

        if self.check[0].get():
            self.root.after(100, self._callback)

    def _disp(self):
        if self.check[1].get():
            self.entries[9]["state"] = "normal"
        else:
            self.entries[9]["state"] = "disabled"

    def _capture(self):
        x1 = self.root.winfo_rootx()
        x2 = x1 + self.root.winfo_reqwidth()
        y1 = self.root.winfo_rooty()
        y2 = y1 + self.root.winfo_reqheight()

        im = grab(bbox=(x1, y1, x2, y2))
        im.show()

    def _plot(self):
        # Parameters
        style.use('seaborn-whitegrid')
        fig1, ax = plt.subplots()
        fig1.canvas.set_window_title('Plot')
        plt.title("Depth Error Chart")
        plt.xlabel("Depth (m)")
        plt.ylabel("Depth error (cm)")
        x1 = []
        y1 = []

        # Light theme
        fig1.patch.set_facecolor('#f5f6f7')  # 212121 Dark
        ax.patch.set_facecolor('#f5f6f7')
        ax.spines['bottom'].set_color('#5c616c')  # FAFAFA Dark
        ax.spines['top'].set_color('#5c616c')
        ax.spines['right'].set_color('#5c616c')
        ax.spines['left'].set_color('#5c616c')
        ax.tick_params(axis='x', colors='#5c616c', which='both')
        ax.tick_params(axis='y', colors='#5c616c', which='both')
        ax.yaxis.label.set_color('#5c616c')
        ax.xaxis.label.set_color('#5c616c')
        ax.title.set_color('#5c616c')

        # Plot
        # Min range to max range
        for x in range(1, 1000):
            y = self._depthErrorCalculator(float(x / 10), float(self.results[1].get() / 1000), float(self.results[0].get()), 0.25)
            x1.append(float(x / 10))
            y1.append(y)

        plt.plot(x1, y1, color="#039BE5")
        plt.legend(['Error'], loc='upper left')

        # Show
        plt.show()


if __name__ == "__main__":
    m = StereoVisionCalculator()
    m.mainloop()
