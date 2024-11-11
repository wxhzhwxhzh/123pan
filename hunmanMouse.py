#!/usr/bin/env python
# -*- coding:utf-8 -*-#

import numpy as np
import math
import random


class bezierTrajectory:

    def _bztsg(self, dataTrajectory):
        lengthOfdata = len(dataTrajectory)

        def staer(x):
            t = ((x - dataTrajectory[0][0]) / (dataTrajectory[-1][0] - dataTrajectory[0][0]))
            y = np.array([0, 0], dtype=np.float64)
            for s in range(len(dataTrajectory)):
                y += dataTrajectory[s] * ((math.factorial(lengthOfdata - 1) / (
                            math.factorial(s) * math.factorial(lengthOfdata - 1 - s))) * math.pow(t, s) * math.pow(
                    (1 - t), lengthOfdata - 1 - s))
            return y[1]

        return staer

    def _type(self, type, x, numberList):
        numberListre = []
        pin = (x[1] - x[0]) / numberList
        if type == 0:
            for i in range(numberList):
                numberListre.append(i * pin)
            if pin >= 0:
                numberListre = numberListre[::-1]
        elif type == 1:
            for i in range(numberList):
                numberListre.append(1 * ((i * pin) ** 2))
            numberListre = numberListre[::-1]
        elif type == 2:
            for i in range(numberList):
                numberListre.append(1 * ((i * pin - x[1]) ** 2))

        elif type == 3:
            dataTrajectory = [np.array([0,0]), np.array([(x[1]-x[0])*0.8, (x[1]-x[0])*0.6]), np.array([x[1]-x[0], 0])]
            fun = self._bztsg(dataTrajectory)
            numberListre = [0]
            for i in range(1,numberList):
                numberListre.append(fun(i * pin) + numberListre[-1])
            if pin >= 0:
                numberListre = numberListre[::-1]
        numberListre = np.abs(np.array(numberListre) - max(numberListre))
        biaoNumberList = ((numberListre - numberListre[numberListre.argmin()]) / (
                    numberListre[numberListre.argmax()] - numberListre[numberListre.argmin()])) * (x[1] - x[0]) + x[0]
        biaoNumberList[0] = x[0]
        biaoNumberList[-1] = x[1]
        return biaoNumberList

    def getFun(self, s):
        '''

        :param s: 传入P点
        :return: 返回公式
        '''
        dataTrajectory = []
        for i in s:
            dataTrajectory.append(np.array(i))
        return self._bztsg(dataTrajectory)

    def simulation(self, start, end, le=1, deviation=0, bias=0.5):
        '''

        :param start:开始点的坐标 如 start = [0, 0]
        :param end:结束点的坐标 如 end = [100, 100]
        :param le:几阶贝塞尔曲线，越大越复杂 如 le = 4
        :param deviation:轨迹上下波动的范围 如 deviation = 10
        :param bias:波动范围的分布位置 如 bias = 0.5
        :return:返回一个字典equation对应该曲线的方程，P对应贝塞尔曲线的影响点
        '''
        start = np.array(start)
        end = np.array(end)
        cbb = []
        if le != 1:
            e = (1 - bias) / (le - 1)
            cbb = [[bias + e * i, bias + e * (i + 1)] for i in range(le - 1)]

        dataTrajectoryList = [start]

        t = random.choice([-1, 1])
        w = 0
        for i in cbb:
            px1 = start[0] + (end[0] - start[0]) * (random.random() * (i[1] - i[0]) + (i[0]))
            p = np.array([px1, self._bztsg([start, end])(px1) + t * deviation])
            dataTrajectoryList.append(p)
            w += 1
            if w >= 2:
                w = 0
                t = -1 * t

        dataTrajectoryList.append(end)
        return {"equation": self._bztsg(dataTrajectoryList), "P": np.array(dataTrajectoryList)}

    def trackArray(self, start, end, numberList, le=1, deviation=10, bias=0.5, type=3, cbb=2, yhh=10):
        '''

        :param start:开始点的坐标 如 start = [0, 0]
        :param end:结束点的坐标 如 end = [100, 100]
        :param numberList:返回的数组的轨迹点的数量 numberList = 150
        :param le:几阶贝塞尔曲线，越大越复杂 如 le = 4
        :param deviation:轨迹上下波动的范围 如 deviation = 10
        :param bias:波动范围的分布位置 如 bias = 0.5
        :param type:0表示均速滑动，1表示先慢后快，2表示先快后慢，3表示先慢中间快后慢 如 type = 1
        :param cbb:在终点来回摆动的次数
        :param yhh:在终点来回摆动的范围
        :return:返回一个字典trackArray对应轨迹数组，P对应贝塞尔曲线的影响点
        '''
        s = []
        fun = self.simulation(start, end, le, deviation, bias)
        w = fun['P']
        fun = fun["equation"]
        if cbb != 0:
            numberListOfcbb = round(numberList*0.2/(cbb+1))
            numberList -= (numberListOfcbb*(cbb+1))

            xTrackArray = self._type(type, [start[0], end[0]], numberList)
            for i in xTrackArray:
                s.append([i, fun(i)])
            dq = yhh/cbb
            kg = 0
            ends = np.copy(end)
            for i in range(cbb):
                if kg == 0:
                    d = np.array([end[0] + (yhh-dq*i), ((end[1]-start[1])/(end[0]-start[0]))*(end[0]+(yhh-dq*i)) +(end[1]-((end[1]-start[1])/(end[0]-start[0]))*end[0])   ])
                    kg = 1
                else:
                    d = np.array([end[0] - (yhh - dq * i), ((end[1]-start[1])/(end[0]-start[0]))*(end[0]-(yhh-dq*i)) +(end[1]-((end[1]-start[1])/(end[0]-start[0]))*end[0])  ])
                    kg = 0
                print(d)
                y = self.trackArray(ends, d, numberListOfcbb, le=2, deviation=0, bias=0.5, type=0, cbb=0, yhh=10)
                s += list(y['trackArray'])
                ends = d
            y = self.trackArray(ends, end, numberListOfcbb, le=2, deviation=0, bias=0.5, type=0, cbb=0, yhh=10)
            s += list(y['trackArray'])

        else:
            xTrackArray = self._type(type, [start[0], end[0]], numberList)
            for i in xTrackArray:
                s.append([i, fun(i)])
        return {"trackArray": np.array(s), "P": w}

class humanMouse:
    def __init__(self):
        self.aa = bezierTrajectory()

    def getRandomTrackArray(self, 水平移动,垂直移动=10,分割=10):
        randomTrackArray = self.aa.trackArray(start=[0,0],end=[水平移动,垂直移动],numberList =分割,deviation = 5,bias = 0.5,type = 3,le=3)
        return randomTrackArray['trackArray']
    def getRandomTrackSpacingArray(self, 水平移动:int=100,垂直移动:int=10,分割:int=10):
        """
        生成一个随机轨迹间距数组。

        参数:
        - 水平移动 (int): 水平方向的移动距离。
        - 垂直移动 (int, 可选): 垂直方向的移动距离，默认为10。
        - 分割 (int, 可选): 分割的段数，即生成的数组长度，默认为10。

        返回:
        - list: 包含随机轨迹间距的数组。
        """
        # 生成一个随机轨迹数组
        randomTrackArray = self.aa.trackArray(start=[0,0],end=[水平移动,垂直移动],numberList =分割,deviation = 5,bias = 0.5,type = 3,le=3)
        # 获取轨迹数组
        bb= randomTrackArray['trackArray']
        # 将第一个点的水平和垂直间距设置为0
        cc=[[0,0]]
        # 获取数组长度
        length = len(bb)
        # 遍历数组，计算每个点与前一个点的水平和垂直间距
        for i in range(1,length):
            x1,y1=bb[i-1]
            x2,y2=bb[i]
            cc.append([float(x2-x1),float(y2-y1)])

        return cc







if __name__ == '__main__':
    # aa=bezierTrajectory()
    # 坐标数组=aa.trackArray(start=[0,10],end=[50,10],numberList =30,deviation = 5,bias = 0.5,type = 3,le=3)
    # print('生产的轨迹坐标',坐标数组['trackArray'])
    from pprint import pprint



    # zuobiao=坐标数组['trackArray']
    轨迹点坐标数组=humanMouse().getRandomTrackArray(147,10,20)
    pprint(轨迹点坐标数组)



    import matplotlib.pyplot as plt

    # 假设这是你的坐标列表
    coordinates = 轨迹点坐标数组

    # 分别提取 x 和 y 坐标
    x, y = zip(*coordinates)

    # 创建一个散点图
    plt.scatter(x, y)

    # 设置图表标题和标签
    plt.title('Scatter Plot of Coordinates')
    plt.xlabel('X axis')
    plt.ylabel('Y axis')

    # 显示网格
    plt.grid()

    # 展示图形
    plt.show()

