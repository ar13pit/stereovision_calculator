from __future__ import (print_function, division, absolute_import,
                        unicode_literals)
from builtins import *
import math
import tkinter as tk
import matplotlib.pyplot as plt


class StereoVisionCalculator(object):
    def __init__(self):
        self.root = tk.Tk()
        self.root.config(bg='#212121')
        self.root.tk_setPalette(background='#212121',
                                activebackground='#212121',
                                fg='#fff',
                                activeforeground='#ccc')

        self.root.title("StereoVision Calculator")
        self.var1 = tk.IntVar(self.root)
        self.pmenu1 = tk.StringVar(self.root)
        self.pmenu2 = tk.StringVar(self.root)
        self.l1 = tk.DoubleVar(self.root)
        self.l2 = tk.DoubleVar(self.root)
        self.l3 = tk.DoubleVar(self.root)
        self.l4 = tk.DoubleVar(self.root)

        # Entries
        self.entries = [tk.Entry(self.root) for _ in range(8)]

        for row, entry in enumerate(self.entries, start=1):
            entry.grid(row=row, column=1, sticky="W")

        # Text entries
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
                       variable=self.var1,
                       command=self._callback,
                       selectcolor='#212121').grid(row=0, sticky="W")
        tk.Button(self.root, text="Help", width=10,
                  command=self._help).grid(row=0, column=2, sticky="W")
        tk.Button(self.root,
                  text="Calculate",
                  width=10,
                  command=self._callback).grid(row=13, sticky="W")
        tk.Button(self.root, text="Plot", width=10,
                  command=self._plot).grid(row=13, column=2, sticky="W")

        # Dropdown menus
        choices1 = {'mm', 'inch'}
        self.pmenu1.set('mm')  # set the default option
        self.popupMenu1 = tk.OptionMenu(self.root, self.pmenu1, *choices1)
        self.popupMenu1.grid(row=1, column=2, sticky="W")

        choices2 = {'Horizontal', 'Vertical', 'Diagonal'}
        self.pmenu2.set('Horizontal')  # set the default option
        self.popupMenu2 = tk.OptionMenu(self.root, self.pmenu2, *choices2)
        self.popupMenu2.grid(row=4, column=2, sticky="W")

        self.popupMenu1.config(highlightthickness='0')
        self.popupMenu2.config(highlightthickness='0')

        # Result entries
        tk.Label(self.root, textvariable=self.l1).grid(row=9,
                                                       column=1,
                                                       sticky="W")
        tk.Label(self.root, textvariable=self.l2).grid(row=10,
                                                       column=1,
                                                       sticky="W")
        tk.Label(self.root, textvariable=self.l3).grid(row=11,
                                                       column=1,
                                                       sticky="W")
        tk.Label(self.root, textvariable=self.l4).grid(row=12,
                                                       column=1,
                                                       sticky="W")

    def mainloop(self):
        self.root.mainloop()

    def _focalLengthCalculator(self, sensor_size, img_width, img_height,
                               focal_fov):
        """
        Function to calculate the focal length of the imaging sensor given:
        :param: sensor_size: The diagonal sensor size in inch or mm, eg: 2/3in
        :param: img_width: The amount of pixels in width
        :param: img_height:  The amount of pixels in height
        :param: focal_fov: The horizontal/vertical/diagonal field of view of the lens, eg: 70h
        """
        if ('inch') in sensor_size:
            sensor_size = sensor_size.replace('inch', '')
            sensor_size = float(sensor_size) * 25.4
        else:
            sensor_size = float(sensor_size.replace('mm', ''))

        ratio = img_height / img_width
        sensor_width_mm = math.sqrt(sensor_size * sensor_size /
                                    (1.0 + ratio * ratio))
        sensor_height_mm = math.sqrt(sensor_size * sensor_size /
                                     (1.0 + 1.0 / (ratio * ratio)))

        roi_width_mm = sensor_width_mm  # * roi_width / img_width
        roi_height_mm = sensor_height_mm  # * roi_height / img_height
        roi_diagonal_mm = math.sqrt(roi_height_mm * roi_height_mm +
                                    roi_width_mm * roi_width_mm)

        fov = float(focal_fov[:-1]) / 180 * math.pi
        atanInner = math.tan(fov * 0.5)

        if ('H') in focal_fov:
            f_mm = roi_width_mm / (2 * atanInner)
        elif ('V') in focal_fov:
            f_mm = roi_height_mm / (2 * atanInner)
        elif ('D') in focal_fov:
            f_mm = roi_diagonal_mm / (2 * atanInner)

        pixel_size_mm = roi_width_mm / img_width

        f_pixel = f_mm / pixel_size_mm

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

        if (self.entries[0].get() and self.entries[1].get()
                and self.entries[2].get() and self.entries[3].get()
                and self.pmenu1.get() and self.pmenu2.get()):
            sensor_size = self.entries[0].get() + self.pmenu1.get()
            img_width = int(self.entries[1].get())
            img_height = int(self.entries[2].get())
            focal_fov = self.entries[3].get() + self.pmenu2.get()[0]

            f_mm, f_pixel = self._focalLengthCalculator(
                sensor_size, img_width, img_height, focal_fov)
            self.l1.set(round(f_mm, 2))

            if (self.entries[4].get() and self.entries[5].get()
                    and self.entries[6].get() and self.entries[7].get()):
                focal_length = f_pixel
                max_depth = float(self.entries[4].get())
                max_depth_error = float(self.entries[5].get())
                min_disparity_measured = int(self.entries[6].get())
                disparity_range = int(self.entries[7].get())

                baseline, min_depth, max_disparity_error = self._baselineCalculator(
                    min_disparity_measured, max_depth, max_depth_error,
                    disparity_range, focal_length)

                self.l2.set(round(baseline * 1000, 2))
                self.l3.set(round(min_depth * 100, 2))
                self.l4.set(round(max_disparity_error, 4))

        if self.var1.get():
            self.root.after(100, self._callback)

    def _help(self):
        m1 = tk.Tk()
        m1.mainloop()

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
        x = self.l3.get()
        y = self.l4.get()
        plt.plot(x, y)

        # Show
        plt.show()


if __name__ == "__main__":
    m = StereoVisionCalculator()
    m.mainloop()
