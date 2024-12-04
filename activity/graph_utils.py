import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import numpy as np


class LiveTimeseriesGraph:
    def __init__(self, time_window, dt, title="", y_label="", legend_label=""):
        plt.ion()
        self.time_window = time_window
        self.n_times = (int)(time_window / dt)
        self.time = np.zeros(self.n_times, dtype=np.float32)
        self.vals = np.zeros(self.n_times, dtype=np.float32)
        self.dt = dt

        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        (self.p1,) = self.ax.plot([], [], label=legend_label)
        self.ax.set_xlabel("Time (ps)")
        self.ax.set_ylabel(y_label)
        self.ax.set_xlim(0, time_window)

        self.ax.set_title(title)
        self.ax.legend()
        self.ax.grid(True)

        self.i = 0

    def update(self, time, val):
        self.idx = self.i % self.n_times
        self.time[self.idx] = time
        self.vals[self.idx] = val

        if self.idx == 0:
            self.ax.set_xlim(
                self.time[self.idx],
                self.time[self.idx] + self.n_times * self.dt,
            )

        valid_indices = range(0, self.idx + 1)

        self.ax.set_ylim(
            np.min(self.vals[valid_indices]) - 1, np.max(self.vals) + 1
        )

        self.p1.set_data(self.time[valid_indices], self.vals[valid_indices])
        clear_output(wait=True)
        display(self.fig)

        self.i += 1
