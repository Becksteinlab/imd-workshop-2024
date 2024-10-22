import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
plt.ion()

class dynamicPlot:
    def __init__(self,xRange,xLabel,yRange1,yRange2,yLabel,titles):
        self.xRange = xRange
        self.xLabel = xLabel
        self.yRange1 = yRange1
        self.yRange2 = yRange2
        self.yLabel = yLabel
        self.titles = titles
        self.prepare()

    def prepare(self):
        #Set up plot
        self.figure, (self.ax1, self.ax2) = plt.subplots(2,1)
        plt.subplots_adjust(hspace = 0.75)
        self.l1, = self.ax1.plot([],[], color = "red", linewidth = 2.0)
        self.ax1.set_xlim(self.xRange[0], self.xRange[1])
        self.ax1.set_ylim(self.yRange1[0], self.yRange1[1])
        self.ax1.set_xlabel(self.xLabel)
        self.ax1.set_ylabel(self.yLabel)
        self.ax1.set_title(self.titles[0])
        self.l2, = self.ax2.plot([],[], color = "blue", linewidth = 2.0)
        self.ax2.set_xlim(self.xRange[0], self.xRange[1])
        self.ax2.set_ylim(self.yRange2[0], self.yRange2[1])
        self.ax2.set_xlabel(self.xLabel)
        self.ax2.set_ylabel(self.yLabel)
        self.ax2.set_title(self.titles[1])

    def update(self, xData, yData1, yData2):
        #Update data (with the new _and_ the old points)
        self.l1.set_xdata(xData)
        self.l1.set_ydata(yData1)
        self.l2.set_xdata(xData)
        self.l2.set_ydata(yData2)
        #We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
    
    def update2(self, xData, yData1, yData2, xRange):
        #Update x-axis limits
        self.xRange = xRange
        self.ax1.set_xlim(self.xRange[0], self.xRange[1])
        self.ax2.set_xlim(self.xRange[0], self.xRange[1])
        #Update data (with the new _and_ the old points)
        self.l1.set_xdata(xData)
        self.l1.set_ydata(yData1)
        self.l2.set_xdata(xData)
        self.l2.set_ydata(yData2)
        #We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
