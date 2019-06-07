from __future__ import (print_function, division, absolute_import,
                        unicode_literals)
from builtins import *
import math
import tkinter as tk
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
        self.root = tk.Tk()
        self.root.tk_setPalette(background='#212121',
                                activebackground='#212121',
                                fg='#fff',
                                activeforeground='#ccc')
        self.root.title("StereoVision Calculator")

        self.var = [tk.IntVar(self.root) for _ in range(2)]  # Checkboxes
        self.pmenu = [tk.StringVar(self.root) for _ in range(2)]  # Popup menus
        self.l = [tk.DoubleVar(self.root) for _ in range(4)]  # Results
        self.entries = [tk.Entry(self.root) for _ in range(8)]  # Entries
        self.entries[7]["state"] = "disabled"

        # Entries
        for row, entry in enumerate(self.entries, start=1):
            entry.grid(row=row, column=1, sticky="W")

        # Text
        textlist = [
            'Sensor size', 'Resolution width', 'Resolution height',
            'Focal FoV', 'Max depth', 'Max depth error', 'Min disparity',
            'Disparity range', 'Focal length', 'Baseline',
            'Minimum measurable depth',
            'Max allowed disparity error \n from calibration'
        ]
        symbollist = [
            'pixels', 'pixels', '', 'm', 'm', 'pixels', 'pixels', 'mm', 'mm',
            'cm', 'pixels'
        ]

        i = 0
        while i < len(textlist):
            tk.Label(self.root, text=textlist[i]).grid(row=i + 1, sticky="W")
            i += 1

        j = 0
        while j < len(symbollist):
            tk.Label(self.root, text=symbollist[j]).grid(row=j + 2,
                                                         column=2,
                                                         sticky="W")
            j += 1

        # Buttons
        tk.Checkbutton(self.root,
                       text="Auto calculate",
                       variable=self.var[0],
                       command=self._callback,
                       selectcolor='#212121').grid(row=0, sticky="W")
        tk.Checkbutton(self.root,
                       text="",
                       variable=self.var[1],
                       command=self._disp,
                       selectcolor='#212121').grid(row=8, column=1, sticky="E")
        tk.Button(self.root,
                  text="Capture",
                  width=12,
                  command=self._capture).grid(row=0, column=2, sticky="W")
        tk.Button(self.root,
                  text="Calculate",
                  width=12,
                  command=self._callback).grid(row=13, sticky="W")
        tk.Button(self.root,
                  text="Plot",
                  width=12,
                  command=self._plot).grid(row=13, column=2, sticky="W")

        # Dropdown menus
        choices1 = {'mm', 'inch'}
        self.pmenu[0].set('mm')  # set the default option
        self.popupMenu0 = tk.OptionMenu(self.root, self.pmenu[0], *choices1)
        self.popupMenu0.grid(row=1, column=2, sticky="W")

        choices2 = {'Horizontal', 'Vertical', 'Diagonal'}
        self.pmenu[1].set('Horizontal')  # set the default option
        self.popupMenu1 = tk.OptionMenu(self.root, self.pmenu[1], *choices2)
        self.popupMenu1.grid(row=4, column=2, sticky="W")

        self.popupMenu0.config(highlightthickness='0', width=12)
        self.popupMenu1.config(highlightthickness='0', width=12)

        # Results
        tk.Label(self.root, textvariable=self.l[0]).grid(row=9,
                                                         column=1,
                                                         sticky="W")
        tk.Label(self.root, textvariable=self.l[1]).grid(row=10,
                                                         column=1,
                                                         sticky="W")
        tk.Label(self.root, textvariable=self.l[2]).grid(row=11,
                                                         column=1,
                                                         sticky="W")
        tk.Label(self.root, textvariable=self.l[3]).grid(row=12,
                                                         column=1,
                                                         sticky="W")

    def mainloop(self):
        self.root.mainloop()

    def _focalLengthCalculator(self, sensor_size, size, img_width, img_height,
                               focal_fov, fov_type):
        """
        Function to calculate the focal length of the imaging sensor given:
        :param: sensor_size: The diagonal sensor size in inch or mm, eg: 2/3in
        :param: img_width: The amount of pixels in width
        :param: img_height:  The amount of pixels in height
        :param: focal_fov: The horizontal/vertical/diagonal field of view of the lens, eg: 70h
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
        return depth_measured - depth

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
            self.l[0].set(round(f_mm, 2))

            if (self.entries[4].get() and self.entries[5].get() and
                    self.entries[6].get() and self.entries[7].get()):
                focal_length = f_pixel
                max_depth = float(self.entries[4].get())
                max_depth_error = float(self.entries[5].get())
                min_disparity_measured = int(self.entries[6].get())
                disparity_range = int(self.entries[7].get())

                baseline, min_depth, max_disparity_error = self._baselineCalculator(
                    min_disparity_measured, max_depth, max_depth_error,
                    disparity_range, focal_length)

                self.l[1].set(round(baseline * 1000, 2))
                self.l[2].set(round(min_depth * 100, 2))
                self.l[3].set(round(max_disparity_error, 4))

        if self.var[0].get():
            self.root.after(100, self._callback)

    def _disp(self):
        if self.var[1].get():
            self.entries[7]["state"] = "normal"
        else:
            self.entries[7]["state"] = "disabled"

    def _capture(self):
        x1 = self.root.winfo_rootx()
        x2 = x1 + self.root.winfo_reqwidth()
        y1 = self.root.winfo_rooty()
        y2 = y1 + self.root.winfo_reqheight()

        im = grab(bbox=(x1, y1, x2, y2))
        im.show()

    def _plot(self):
        # Parameters
        fig1, ax = plt.subplots()
        fig1.canvas.set_window_title('Plot')
        plt.title("Depth Error Chart")
        plt.xlabel("Depth (m)")
        plt.ylabel("Depth error (cm)")

        # Dark theme
        fig1.patch.set_facecolor('#212121')
        ax.patch.set_facecolor('#212121')
        ax.spines['bottom'].set_color('#FAFAFA')
        ax.spines['top'].set_color('#FAFAFA')
        ax.spines['right'].set_color('#FAFAFA')
        ax.spines['left'].set_color('#FAFAFA')
        ax.tick_params(axis='x', colors='#FAFAFA', which='both')
        ax.tick_params(axis='y', colors='#FAFAFA', which='both')
        ax.yaxis.label.set_color('#FAFAFA')
        ax.xaxis.label.set_color('#FAFAFA')
        ax.title.set_color('#FAFAFA')

        # Plot
        x = self.l[2].get()
        y = self.l[3].get()
        plt.plot(x, y)

        # Show
        plt.show()


if __name__ == "__main__":
    m = StereoVisionCalculator()
    m.mainloop()
