from __future__ import (print_function, division, absolute_import,
                        unicode_literals)
from builtins import *
import math
import tkinter as tk
from tkinter import ttk
from ttkthemes import ThemedTk
import matplotlib.pyplot as plt
from pyscreenshot import grab


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

        # Baseline calculator parameters
        self.max_depth = None
        self.max_depth_error = None
        self.max_disparity = None
        self.calibration_disparity_error = None
        self.min_disparity = None
        self.baseline = None
        self.min_depth = None

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
        self.pmenu = [tk.StringVar(self.root) for _ in range(2)]  # Popup menus
        self.results = [tk.DoubleVar(self.root) for _ in range(5)]  # Results
        self.entries = [ttk.Entry(self.root) for _ in range(12)]  # Entries
        self.entries[9]["state"] = "disabled"
        self.entries[9]["width"] = 14

        self.results[3].set('0 - 0')
        self.results[4].set('0')

        # Entries
        for row, entry in enumerate(self.entries, start=1):
            entry.grid(row=row, column=1, sticky="W")

        textlist = {
            'Sensor size': '',
            'Resolution width': 'pixels',
            'Resolution height': 'pixels',
            'Focal FoV': '',
            'Min baseline': 'mm',
            'Max baseline': 'mm',
            'Min depth': 'mm',
            'Performance depth': 'm',
            'Performance depth error': 'cm',
            'Performance disparity': 'pixels',
            'Max disparity': 'pixels',
            'Calibration disparity error': 'pixel',
            'Focal length': 'mm',
            'Baseline': 'mm',
            'Max depth': 'm',
            'Depth resolution': 'pixels',
            'Depth FoV': 'deg'
        }

        for x, (val, key) in enumerate(textlist.items()):
            ttk.Label(self.root, text=val).grid(
                row=x + 1, sticky="W", pady=5)
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
                   command=self._callback).grid(row=18, sticky="W")
        ttk.Button(self.root,
                   text="Plot",
                   width=12,
                   command=self._plot).grid(row=18, column=2, sticky="W")

        # Dropdown menus
        choices1 = {'', 'mm', 'inch'}
        self.pmenu[0].set('mm')  # set the default option
        self.popupMenu0 = ttk.OptionMenu(self.root, self.pmenu[0], *choices1)
        self.popupMenu0.grid(row=1, column=2, sticky="W")

        choices2 = {'', 'Horizontal', 'Vertical', 'Diagonal'}
        self.pmenu[1].set('Horizontal')  # set the default option
        self.popupMenu1 = ttk.OptionMenu(self.root, self.pmenu[1], *choices2)
        self.popupMenu1.grid(row=4, column=2, sticky="W")

        # Results
        for x in range(0, len(self.results)):
            ttk.Label(self.root, textvariable=self.results[x]).grid(
                row=x + 13, column=1, sticky="W")

        col_count, row_count = self.root.grid_size()

        for col in range(col_count):
            self.root.grid_columnconfigure(col, pad=2)

    def mainloop(self):
        self.root.mainloop()

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

        return f_mm, f_pixel

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

    def _callback(self):

        if (self.entries[0].get() and self.entries[1].get() and
                self.entries[2].get() and self.entries[3].get() and
                self.pmenu[0].get() and self.pmenu[1].get()):
            sensor_size = float(self.entries[0].get())
            size = self.pmenu[0].get()
            img_width = int(self.entries[1].get())
            img_height = int(self.entries[2].get())
            focal_fov = float(self.entries[3].get())
            fov_type = self.pmenu[1].get()

            f_mm, f_pixel = self._focalLengthCalculator(
                sensor_size, size, img_width, img_height, focal_fov, fov_type)
            self.results[0].set(round(f_mm, 2))

            if (self.entries[4].get() and self.entries[5].get() and
                    self.entries[6].get() and self.entries[7].get()):
                focal_length = f_pixel
                max_depth = float(self.entries[4].get())
                max_depth_error = float(self.entries[5].get())
                calibration_disparity_error = float(self.entries[6].get())
                max_disparity = int(self.entries[7].get())
                min_disparity = 1 if not self.entries[8].get() else int(
                    self.entries[8].get())

                baseline, min_depth = self._baselineCalculator(
                    focal_length, max_depth, max_depth_error,
                    max_disparity, calibration_disparity_error,
                    min_disparity)

                self.results[1].set(round(baseline * 1000, 2))
                self.results[2].set(round(min_depth * 100, 2))

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

    def _xfrange(self, start, stop, step):
        i = 0
        while start + i * step < stop:
            yield start + i * step
            i += 1

    def _plot(self):
        # Parameters
        fig1, ax = plt.subplots()
        fig1.canvas.set_window_title('Plot')
        plt.title("Depth Error Chart")
        plt.xlabel("Depth (m)")
        plt.ylabel("Depth error (cm)")
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 10])
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
        for x in self._xfrange(1.0, 100.0, 0.1):
            y = self._depthErrorCalculator(x, float(self.results[1].get() / 1000), float(self.results[0].get()), 0.25)
            x1.append(x)
            y1.append(y)

        print("X: " + str(x1) + "\tY: " + str(y1))

        plt.plot(x1, y1, color="b")
        plt.legend(('Error'), loc='upper left')

        # Show
        plt.show()


if __name__ == "__main__":
    m = StereoVisionCalculator()
    m.mainloop()
