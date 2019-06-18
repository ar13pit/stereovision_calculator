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

    Row = namedtuple('Row', ['name', 'io', 'activate', 'units'])

    class RowElement(object):
        def __init__(self):
            self.name = None
            self.io = None
            self.activate = None
            self.units = None

        def setrow(self, row):
            self.name.grid(row=row)
            self.io.grid(row=row)
            self.activate.grid(row=row)
            self.units.grid(row=row)

    def _rowPropertiesToGUI(self, master, row_prop):
        """
        Method to convert row_prop of type Row() into a GUI element
        :param: master (The master of tk element)
        :param: row_prop (An instance of Row())
        :return: An instance of RowElement
        """
        row = RowElement()

        # Create the name label
        row.name = ttk.Label(master, text=row_prop.name)
        row.name.grid(column=0, sticky="W", pady=5)

        # Create io Entry var
        if row_prop.io:
            row.io = ttk.Entry(master)
        else:
            row.io = tk.StringVar(master)
        row.io.grid(column=1, sticky="W")

        # Create units Label/OptionMenu var
        if row_prop.units:
            if isinstance(row_prop.units, list):
                row.units = ttk.OptionMenu(master, tk.StringVar(master,
                    row_prop.units[0]), *row_prop.units)
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

        sensor_size_units = ['mm', 'in']
        fov_type = ['Horizontal', 'Vertical', 'Diagonal']
        depth_units = ['mm', 'cm', 'm', 'in', 'ft']
        row_properties = [
            Row('Sensor size', True, False, sensor_size_units),
            Row('Resolution width', True, False, 'px'),
            Row('Resolution height', True, False, 'px'),
            Row('Focal FoV', True, False, fov_type),
            Row('Performance depth', True, False, depth_units),
            Row('Performance depth error', True, False, depth_units),
            Row('Performance disparity', True, True, 'px'),
            Row('Max disparity', True, False, 'px'),
            Row('Calibration disparity error', True, False, 'px'),
            Row('Focal length', False, False, 'mm'),
            Row('Baseline', False, False, 'cm'),
            Row('Min depth', False, False, 'm'),
            Row('Max depth', False, False, 'm'),
            Row('Depth resolution', False, False, 'px'),
            Row('Depth FoV', False, False, 'deg')
        ]

        self.check = [tk.IntVar(self.root) for _ in range(2)]  # Checkboxes
        self.pmenu = [tk.StringVar(self.root) for _ in range(4)]  # Popup menus
        self.results = [tk.DoubleVar(self.root) for _ in range(6)]  # Results
        self.entries = [ttk.Entry(self.root) for _ in range(9)]  # Entries
        self.entries[6]["state"] = "disabled"
        self.entries[6]["width"] = 14

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
            l = ttk.Label()
            l.master = self.root
            l["text"] = val
            l.grid(row=x+1, sticky="W", pady=5)
            #ttk.Label(self.root, text=val).grid(
            #    row=x + 1, sticky="W", pady=5)
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
                        command=self._disp).grid(row=7, column=1, sticky="E")
        ttk.Button(self.root,
                   text="Capture",
                   width=12,
                   command=self._capture).grid(row=0, column=2, sticky="W")
        ttk.Button(self.root,
                   text="Calculate",
                   width=12,
                   command=self._callback).grid(row=16, sticky="W")
        ttk.Button(self.root,
                   text="Plot",
                   width=12,
                   command=self._plot).grid(row=16, column=2, sticky="W")

        # Dropdown menus
        choices1 = {'', 'mm', 'in'}
        choices2 = {'', 'Horizontal', 'Vertical', 'Diagonal'}
        choices3 = {'', 'mm', 'cm', 'm', 'in', 'ft'}
        self.popupMenu = []
        for x in range(0, 4):
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
                row=x + 10, column=1, sticky="W")

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
        sensor_width_mm = math.sqrt(sensor_size * sensor_size /
                                    (1.0 + ratio * ratio))
        sensor_height_mm = math.sqrt(sensor_size * sensor_size /
                                     (1.0 + 1.0 / (ratio * ratio)))

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

            pixel_size_mm = roi_width_mm / img_width
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
        self._min_depth = disparityToDepth(self._baseline, self._focal_length,
                self.perf_disp_max)

    def _depthErrorCalculator(self, depth):
        """
        Method to calculate the max_depth_error for a given depth for
        pre-determined baseline and focal length
        :param: depth
        :return: max_depth_error
        """
        disparity_real = depthToDisparity(self._baseline, self._focal_length,
                depth)
        disparity_measured = disparity_real + self.perf_disp_calibration_error

        depth_measured = disparityToDepth(self._baseline, self._focal_length,
                disparity_measured)
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

            self.sensor_size = float(self.entries[0].get())
            size = self.pmenu[0].get()
            self.img_width = int(self.entries[1].get())
            self.img_height = int(self.entries[2].get())
            self.focal_fov = float(self.entries[3].get())
            fov_type = self.pmenu[1].get()

            f_mm, roi_width_mm, roi_height_mm = self._focalLengthCalculator(
                size, fov_type)
            self.results[0].set(round(f_mm, 2))

            if (self.entries[4].get() and self.entries[5].get() and
                    self.entries[6].get() and self.entries[7].get() and
                    self.entries[10].get()):
                self.perf_depth = calculateToMeter(
                    float(self.entries[4].get()), self.pmenu[7].get)
                self.perf_depth_error = calculateToMeter(
                    float(self.entries[5].get()), self.pmenu[8].get)
                self.perf_disp = 1 if not self.entries[9].get() else int(
                    self.entries[9].get())
                self.perf_disp_max = int(self.entries[10].get())
                self.perf_disp_calibration_error = float(self.entries[11].get())

                self._baselineCalculator()

                d_max = self.perf_disp_max - 1
                depth_fov, depth_res = self._depthCalculator(
                    self.img_width, self.img_height, roi_width_mm,
                    roi_height_mm, self.img_width, d_max, f_mm)
                self.results[1].set(round(self._baseline * 1000, 2))
                self.results[2].set(round(self._min_depth * 100, 2))
                self.results[4].set(depth_res)
                self.results[5].set(depth_fov)

        if self.check[0].get():
            self.root.after(100, self._callback)

    def _disp(self):
        if self.check[1].get():
            self.entries[6]["state"] = "normal"
        else:
            self.entries[6]["state"] = "disabled"

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
