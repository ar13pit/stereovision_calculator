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
from collections import namedtuple


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
    return (baseline * focal_length / disparity)


def depthToDisparity(baseline, focal_length, depth):
    """
    Function to calculate the disparity between image points (in px) given:
    :param: baseline: The distance between the camera centers (in m)
    :param: focal_length: The focal length of the system (in px)
    :param: depth: Depth of the image point (in m)
    """
    return (baseline * focal_length / depth)


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

        # Design limits
        self.perf_depth = None
        self.perf_depth_error = None
        self.perf_disp_max = None
        self.perf_disp = None
        self.perf_disp_calibration_error = None

        # Baseline calculator results
        self._min_depth = None
        self._baseline = None
        self._focal_length = None

        # Initialize the complete GUI
        self._initializeGUI()

    Row = namedtuple('Row', ['var_name', 'name', 'io', 'activate', 'units'])

    class RowElement(object):
        def __init__(self):
            self.name = None
            self.io = None
            self.activate = None
            self.units = None

        def setrow(self, row):
            if self.name:
                self.name.grid(row=row)
            if self.io:
                self.io.grid(row=row)
            if self.activate:
                self.activate.grid(row=row)
            if self.units:
                self.units.grid(row=row)

    def _rowPropertiesToGUI(self, master, row_prop):
        """
        Method to convert row_prop of type Row() into a GUI element
        :param: master (The master of tk element)
        :param: row_prop (An instance of Row())
        :return: An instance of RowElement
        """
        row = StereoVisionCalculator.RowElement()

        # Create the name label
        row.name = ttk.Label(master, text=row_prop.name)
        row.name.grid(column=0, sticky="W", pady=5)

        # Create io Entry var
        if row_prop.io:
            row.io = ttk.Entry(master)
        else:
            row.io = ttk.Label(master, text="0.0")
        row.io.grid(column=1, sticky="W")

        if row_prop.activate:
            row.activate = ttk.Checkbutton(master)
            row.activate.grid(column=1, sticky="E")

        # Create units Label/OptionMenu var
        if row_prop.units:
            if isinstance(row_prop.units, list):
                row.units = ttk.OptionMenu(master, tk.StringVar(master,
                    row_prop.units[1]), *row_prop.units)
            else:
                row.units = tk.Label(master, text=row_prop.units)

            row.units.grid(column=2, sticky="W")

        return row

    def _initializeGUI(self):
        """
        Method to setup the StereoVision Calculator GUI using tkinter
        """
        self.root = ThemedTk(theme="arc")
        self.root.tk_setPalette(background='#f5f6f7')
        self.root.title("StereoVision Calculator")
        self.root.resizable(0, 0)  # Don't allow resize

        sensor_size_units = ['', 'mm', 'in']
        fov_type = ['', 'Horizontal', 'Vertical', 'Diagonal']
        depth_units = ['', 'mm', 'cm', 'm', 'in', 'ft']
        row_properties = [
            StereoVisionCalculator.Row('ui_sensor_size', 'Sensor size', True, False, sensor_size_units),
            StereoVisionCalculator.Row('ui_img_width', 'Image width', True, False, 'px'),
            StereoVisionCalculator.Row('ui_img_height', 'Image height', True, False, 'px'),
            StereoVisionCalculator.Row('ui_focal_fov', 'Focal FoV', True, False, fov_type),
            StereoVisionCalculator.Row('ui_perf_depth', 'Performance depth', True, False, depth_units),
            StereoVisionCalculator.Row('ui_perf_depth_error', 'Performance depth error', True, False, depth_units),
            StereoVisionCalculator.Row('ui_perf_disp', 'Performance disparity', True, True, 'px'),
            StereoVisionCalculator.Row('ui_disp_max', 'Max disparity', True, False, 'px'),
            StereoVisionCalculator.Row('ui_disp_cal_error', 'Calibration disparity error', True, False, 'px'),
            StereoVisionCalculator.Row('ui_focal_length', 'Focal length', False, False, 'mm'),
            StereoVisionCalculator.Row('ui_baseline', 'Baseline', False, False, 'cm'),
            StereoVisionCalculator.Row('ui_depth_min', 'Min depth', False, False, 'm'),
            StereoVisionCalculator.Row('ui_depth_max', 'Max depth', False, False, 'm'),
            StereoVisionCalculator.Row('ui_depth_res', 'Depth resolution', False, False, 'px'),
            StereoVisionCalculator.Row('ui_depth_fov', 'Depth FoV', False, False, 'deg')
        ]

        for row_num, rp in enumerate(row_properties, start=1):
            row_element = self._rowPropertiesToGUI(self.root, rp)
            row_element.setrow(row_num)
            self.__setattr__(rp.var_name, row_element)

        self.ui_perf_disp.io["state"] = tk.DISABLED
        self.ui_perf_disp.io["width"] = 14
        self.ui_perf_disp.activate["command"] = self._disp

        self.ui_depth_res.io["text"] = "0 x 0"
        self.ui_depth_fov.io["text"] = "0° x 0°"

        # Buttons
        self.ui_auto_calculate = ttk.Checkbutton(self.root,
                                                 text="Auto calculate",
                                                 command=self._callback)
        self.ui_auto_calculate.grid(row=0, sticky="W", pady=5)

        self.ui_capture = ttk.Button(self.root,
                                     text="Capture",
                                     width=12,
                                     command=self._capture)
        self.ui_capture.grid(row=0, column=2, sticky="W")

        self.ui_calculate = ttk.Button(self.root,
                                       text="Calculate",
                                       width=12,
                                       command=self._callback)
        self.ui_calculate.grid(row=16, sticky="W")

        self.ui_plot = ttk.Button(self.root,
                                  text="Plot",
                                  width=12,
                                  command=self._plot)
        self.ui_plot.grid(row=16, column=2, sticky="W")

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

    def _focalLengthCalculator(self, size, fov_type):
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
            self.sensor_size = self.sensor_size * 25.4
        else:
            self.sensor_size = self.sensor_size

        ratio = self.img_height / self.img_width
        sensor_width_mm = math.sqrt(self.sensor_size ** 2 /
                                    (1.0 + ratio ** 2))
        sensor_height_mm = math.sqrt(self.sensor_size ** 2 /
                                     (1.0 + 1.0 / (ratio ** 2)))

        roi_width_mm = sensor_width_mm  # * roi_width / img_width
        roi_height_mm = sensor_height_mm  # * roi_height / img_height
        roi_diagonal_mm = math.sqrt(roi_height_mm * roi_height_mm +
                                    roi_width_mm * roi_width_mm)

        fov = self.focal_fov / 180 * math.pi
        atanInner = math.tan(fov * 0.5)

        try:
            if fov_type == 'Horizontal':
                f_mm = roi_width_mm / (2 * atanInner)
            elif fov_type == 'Vertical':
                f_mm = roi_height_mm / (2 * atanInner)
            elif fov_type == 'Diagonal':
                f_mm = roi_diagonal_mm / (2 * atanInner)

            pixel_size_mm = roi_width_mm / self.img_width
            self._focal_length = f_mm / pixel_size_mm
        except ZeroDivisionError:
            f_mm = 0
            self._focal_length = 0

        return f_mm, roi_width_mm, roi_height_mm

    def _baselineCalculator(self):
        """
        Function to calculate the baseline and min depth of the stereo camera given:
            1. Focal length of the lens
            2. Performace depth
            3. Performance depth error
            4. Disparity at performance depth
            5. Calibration disparity error
        """
        disparity = self.perf_disp + self.perf_disp_calibration_error
        depth = self.perf_depth + self.perf_depth_error

        self._baseline = baselineCalculator(self._focal_length, disparity, depth)
        self._min_depth = disparityToDepth(self._baseline, self._focal_length, self.perf_disp_max)

    def _depthErrorCalculator(self, depth):
        """
        Method to calculate the max_depth_error for a given depth for
        pre-determined baseline and focal length
        :param: depth
        :return: max_depth_error
        """
        disparity_real = depthToDisparity(self._baseline, self._focal_length, depth)
        disparity_measured = disparity_real + self.perf_disp_calibration_error

        depth_measured = disparityToDepth(self._baseline, self._focal_length, disparity_measured)
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

        if (self.ui_sensor_size.io.get() and self.ui_img_width.io.get() and
                self.ui_img_height.io.get() and self.ui_focal_fov.io.get()):

            self.sensor_size = float(self.ui_sensor_size.io.get())
            size = self.ui_sensor_size.units["text"]
            self.img_width = int(self.ui_img_width.io.get())
            self.img_height = int(self.ui_img_height.io.get())
            self.focal_fov = float(self.ui_focal_fov.io.get())
            fov_type = self.ui_focal_fov.units["text"]

            f_mm, roi_width_mm, roi_height_mm = self._focalLengthCalculator(
                size, fov_type)
            self.ui_focal_length.io["text"] = round(f_mm, 2)

            if (self.ui_perf_depth.io.get() and self.ui_perf_depth_error.io.get() and
                    self.ui_disp_max.io.get() and
                    self.ui_disp_cal_error.io.get()):

                self.perf_depth = calculateToMeter(
                    float(self.ui_perf_depth.io.get()),
                    self.ui_perf_depth.units["text"])
                self.perf_depth_error = calculateToMeter(
                    float(self.ui_perf_depth_error.io.get()),
                    self.ui_perf_depth_error.units["text"])
                self.perf_disp = 1 if not self.ui_perf_disp.io.get() else int(
                    self.ui_perf_disp.io.get())
                self.perf_disp_max = int(self.ui_disp_max.io.get())
                self.perf_disp_calibration_error = float(self.ui_disp_cal_error.io.get())

                self._baselineCalculator()

                d_max = self.perf_disp_max - 1
                depth_fov, depth_res = self._depthCalculator(
                    self.img_width, self.img_height, roi_width_mm,
                    roi_height_mm, self.img_width, d_max, f_mm)
                self.ui_baseline.io["text"] = round(self._baseline * 100, 2)
                self.ui_depth_min.io["text"] = round(self._min_depth * 100, 2)
                self.ui_depth_res.io["text"] = depth_res
                self.ui_depth_fov.io["text"] = depth_fov

        if self.ui_auto_calculate.instate(["selected"]):
            self.root.after(100, self._callback)

    def _disp(self):
        if self.ui_perf_disp.activate.instate(["selected"]):
            self.ui_perf_disp.io["state"] = tk.NORMAL
        else:
            self.ui_perf_disp.io["state"] = tk.DISABLED

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
            y = self._depthErrorCalculator(float(x / 10))
            x1.append(float(x / 10))
            y1.append(y)

        plt.plot(x1, y1, color="#039BE5")
        plt.legend(['Error'], loc='upper left')

        # Show
        plt.show()


if __name__ == "__main__":
    m = StereoVisionCalculator()
    m.mainloop()
