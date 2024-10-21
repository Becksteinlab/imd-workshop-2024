import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
plt.ion()

class dynamicPlot:
    def __init__(self,xRange,xLabel,yRange,yLabel,title):
        self.xRange = xRange
        self.xLabel = xLabel
        self.yRange = yRange
        self.yLabel = yLabel
        self.title = title
        self.prepare()

    def prepare(self):
        #Set up plot
        self.figure, self.ax = plt.subplots()
        self.l1, = self.ax.plot([],[], color = "red", linewidth = 2.0, label = "total")
        self.l2, = self.ax.plot([],[], color = "blue", linewidth = 2.0, label = "single molecule")
        self.ax.set_xlim(self.xRange[0], self.xRange[1])
        self.ax.set_ylim(self.yRange[0], self.yRange[1])
        self.ax.set_xlabel(self.xLabel)
        self.ax.set_ylabel(self.yLabel)
        self.ax.set_title(self.title)
        self.ax.legend()

    def update(self, xData, yData1, yData2):
        #Update data (with the new _and_ the old points)
        self.l1.set_xdata(xData)
        self.l1.set_ydata(yData1)
        self.l2.set_xdata(xData)
        self.l2.set_ydata(yData2)
        #We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
