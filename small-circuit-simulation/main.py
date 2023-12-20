# 导入tkinter模块，用于创建图形化界面
# from solver import solve
import numpy as np
import pandas as pd
import os
from cmath import rect

# 导入matplotlib库
import matplotlib.pyplot as plt
# 导入numpy库
import numpy as np

# 定义一个函数，用于绘制余弦函数图像


def plot_cosine(faai, a, f):
    # 创建一个数组，表示x轴的范围，从0到2π，步长为0.01
    x = np.arange(0, 2 * np.pi/f, 1/f/100)
    # 计算y轴的值，根据余弦函数的公式S
    y = a * np.cos(2 * np.pi * f * x + faai)
    # 绘制图像，设置线条1\的颜色和标签
    plt.plot(x, y, color="blue",
             label="U = {:.2f} * cos(2 * pi * {:.2f} * t + {:.2f})".format(a, f, faai))
    # 设置图像的标题，使用节点编号
    plt.title("voltage-time ")
    # 设置x轴和y轴的标签和单位
    plt.xlabel("time (s)")
    plt.ylabel("voltage (V)")
    # 显示图例
    plt.legend()
    # 修改：设置x轴和y轴的范围，根据频率和幅值的变化
    # 设置x轴的范围为一个周期，即0到1/f
    plt.xlim(0, 3 / f)
    # 设置y轴的范围为幅值的正负，即-a到a
    plt.ylim(-2*a, 2*a)
    # 显示图像
    plt.show()

import math
def draw_pic(df,V):
    # 输入参数    
    frequency = float(df['freq'].head(1).values[0])
    N1 = int(input("请输入第一个节点编号："))
    N2 = int(input("请输入第二个节点编号："))
    N1=N1-1
    N2=N2-1
    a1=np.real(V)[N1]
    b1=np.imag(V)[N1]
    a2=np.real(V)[N2]
    b2=np.imag(V)[N2]
    if (a2-a1)==0:
        if b2*b1>=0:  
            faai =0
        else:
            faai=math.pi/2 
    else:         
        faai = math.atan((b2-b1)/(a2-a1))
    a = math.sqrt((b2-b1)**2+(a2-a1)**2)
    f = frequency
    # 调用函数
    plot_cosine(faai, a, f)



def initialize_GI(node_num):
	# 初始化电导矩阵G和电流矩阵I
	G = np.zeros((node_num, node_num), dtype=complex)
	I = np.zeros((node_num, 1), dtype=complex)
	return G, I


def set_GI(df, branch_num, frequency, G, I):
	# 设置角频率
	omega = 2 * np.pi * frequency

	# 构建电导矩阵G和电流矩阵I
	for i in range(branch_num):
		# 对每一条支路进行遍历
		start_node = int(df.loc[i, 'start']) - 1
		end_node = int(df.loc[i, 'end']) - 1
		# 节点编号是从1开始编的，但是numpy数组检索需要从0开始，故-1
		component_type = df.loc[i, 'type']
		component_value = float(df.loc[i, 'value'])
		if component_type == 'res':
			G[start_node, start_node] += 1 / component_value
			G[end_node, end_node] += 1 / component_value
			G[start_node, end_node] -= 1 / component_value
			G[end_node, start_node] -= 1 / component_value
		elif component_type == 'cap':
			addmitance = 1j * omega * component_value
			G[start_node, start_node] += addmitance
			G[end_node, end_node] += addmitance
			G[start_node, end_node] -= addmitance
			G[end_node, start_node] -= addmitance
		elif component_type == 'ind':
			addmitance = 1 / (1j * omega * component_value)
			G[start_node, start_node] += addmitance
			G[end_node, end_node] += addmitance
			G[start_node, end_node] -= addmitance
			G[end_node, start_node] -= addmitance
		elif component_type == 'idc':
			modulus = component_value
			angle = np.deg2rad(float(df.loc[i, 'phi']))
			# 计算复数的指数形式
			current = rect(modulus, angle)

			I[start_node, 0] -= current
			I[end_node, 0] += current
		elif component_type == 'udc':
			modulus = component_value
			angle = np.deg2rad(float(df.loc[i, 'phi']))
			voltage = rect(modulus, angle)

			# 给电压源串联很小的虚拟电阻，然后进行电源等效
			r_small = 1e-12

			G[start_node, start_node] += 1/r_small
			G[end_node, end_node] += 1/r_small
			G[start_node, end_node] -= 1/r_small
			G[end_node, start_node] -= 1/r_small

			# 添加电流源
			current = voltage / r_small
			I[start_node, 0] -= current
			I[end_node, 0] += current

	# 删除电导矩阵G和电流矩阵I的最后一行和最后一列（最后一个节点作为参考节点）
	G = np.delete(G, -1, axis=0)
	G = np.delete(G, -1, axis=1)
	I = np.delete(I, -1, axis=0)
	return G, I

# 主函数


def solve(df):
    # 电源方向设置为从起始节点到终止节点
    # 获取节点数和支路数
    node_num = df['start'].max()
    branch_num = df.shape[0]
    frequency = float(df['freq'].head(1).values[0])
     
    # 初始化电导矩阵G和电流矩阵I
    G, I = initialize_GI(node_num)
    # 设置电导矩阵
    G, I = set_GI(df, branch_num, frequency, G, I)
    print("电导矩阵为：")
    print(G)
    print("电流矩阵为：")
    print(I)
	# 求解节点电压
    V = np.linalg.solve(G, I)

	# 加上参考节点电压
    V = np.append(V, [[0]], axis=0)
    print("节点电压为：")
    #print(V)
    for row in V:
        print("[", end="")
        for element in row:
            print("{:.2f}".format(element), end=" ")
        print("]")
    while True:
        draw_pic(df,V)
        a=int(input("是否结束 1 or 0:"))
        if a==1:
            print("结束啦，谢谢你！\n本程序由来自22级钱院自动化 智电 智能 计试的\n刘睿勋  杨弋可  李言诣  杨锦添  张跃骞  李兆兴  王铭之  吴美怡  梁立丹  林圣翔共同完成")
            return

import pandas as pd
import tkinter as tk
# 导入PIL模块，用于处理图片
from PIL import Image, ImageTk
import numpy as np
import sys
sys.setrecursionlimit(5000)

print("#####################################################")
print("\033[1m                        说明                         \033[0m")
print("1.先用鼠标选择按键并放置所有元器件，单一元器件可以多次选择放置位置，按空格键确认\n")
print("2.放置所有元器件之后，点击绘制导线，点击某一元件红色端点，当终端显示“选择成功”，可以继续点击其他点放置导线，\n最终点击某一元件蓝色端点，当终端显示“放置成功”即完成一条导线的布置\n")
print("3.放置所有导线之后点击查询关联矩阵，即可查询电路的关联矩阵\n")
print("#####################################################")
WIDTH = 15
value1, value2 = 0, 0
line_num = 0
freq = 0


# 创建一个窗口对象

window = tk.Tk()
var = tk.IntVar()  # 创建一个整型变量
# 设置窗口标题
window.title("my_multism")
# 设置窗口大小
window.geometry("1300x800")
# 创建一个画布对象，用于绘制网格点和元器件
canvas = tk.Canvas(window, width=1200, height=1200, bg="white")

# 创建一个字典，用于存储元器件的图片对象
images = {}
# 创建一个列表，用于存储当前选中的元器件的类型和位置
selected = []

real_selected = []


# 用于加载元器件的图片并缩放到合适的大小
def load_image(name):

    image = Image.open(name + ".png")
    # 获取图片的宽度和高度
    width, height = image.size
    # 计算缩放比例，使图片的宽度或高度不超过50像素
    scale = min(50 / width, 50 / height)
    # 按照缩放比例重新设置图片的宽度和高度
    width = int(width * scale)
    height = int(height * scale)
    # 缩放图片
    image = image.resize((width, height))
    # 将图片转换为tkinter可用的格式
    image = ImageTk.PhotoImage(image)
    # 返回图片对象
    return image

# 定义一个函数，用于绘制网格点


def draw_grid():
    # 遍历每一行
    for i in range(100):
        # 遍历每一列
        for j in range(100):
            # 计算网格点的坐标
            x = 50 + j * WIDTH
            y = 50 + i * WIDTH
            # 在画布上绘制一个小圆点，表示网格点
            canvas.create_oval(x - 1, y - 1, x + 1, y + 1, fill="black")

# 定义一个函数，用于绘制元器件


def draw_image(name, row, col):
    # 从字典中获取元器件的图片对象，如果不存在，则调用load_image函数加载图片并存入字典中
    image = images.get(name)
    if image is None:
        image = load_image(name)
        images[name] = image
    # 计算元器件的中心坐标
    x = 50 + col * WIDTH
    y = 50 + row * WIDTH
    # 在画布上绘制元器件的图片，并返回图片的ID，用于后续移动或删除操作
    return canvas.create_image(x, y, image=image)


# (step,old_x,old_y,now_x,now_y)
# 0代表第一个元器件
line_draw_step = np.zeros(shape=(1, 5), dtype=int)
# 定义一个函数，用于处理鼠标左键点击事件
# 表示模式（是否进入绘制导线模式）
mode = 0
# 连线标号

# 全局数据
df = pd.DataFrame(columns=["type", "start", "end", "value", "phi", "freq"])


def check_point():
    # print("in check_point")
    # print(line_draw_step)
    # print(real_selected)
    global line_num
    num1 = 0
    for i in range(np.shape(real_selected)[0]):
        # 直接连线
        if line_draw_step[0][1] == real_selected[i][1] and line_draw_step[0][2] == real_selected[i][2] and line_draw_step[0][3] == real_selected[i][3] and line_draw_step[0][4] == real_selected[i][4]:
            # print(1)
            num1 += 0
        elif line_draw_step[0][3] == real_selected[i][3] and line_draw_step[0][4] == real_selected[i][4]:
            # print(3)
            num1 += 3
            incidence_matrix[line_num][i] = -1
        # 起点
        elif line_draw_step[0][1] == real_selected[i][1] and line_draw_step[0][2] == real_selected[i][2]:
            # print(2)
            incidence_matrix[line_num][i] = 1
            num1 += 2
        # 终点

    return num1


def on_click(event):
    global line_num
    global selected  # 声明全局变量selected
    # 如果当前没有选中任何元器件，则不做任何操作，直接返回
    x = event.x
    y = event.y
    if mode == 1:
        if line_draw_step[0][0] == 0:
            line_draw_step[0][1] = (x-50)//WIDTH
            line_draw_step[0][2] = (y-50)//WIDTH
            # 选点不对
            if check_point() == 0:
                line_draw_step[0][1] = 0
                line_draw_step[0][2] = 0
                return
            # 选中红色起始点
            elif check_point() == 2:
                line_draw_step[0][1] = (x-50)//WIDTH
                line_draw_step[0][2] = (y-50)//WIDTH
                line_draw_step[0][0] = 1
                print("已选中元器件红色起点")

        elif line_draw_step[0][0] == 1:

            # 选中蓝色终止点
            line_draw_step[0][3] = (x-50)//WIDTH
            line_draw_step[0][4] = (y-50)//WIDTH
            if check_point() == 5 or check_point() == 3:
                draw_line()
                line_draw_step[0][0] = 0
                line_draw_step[0][1] = 0
                line_draw_step[0][2] = 0
                line_draw_step[0][3] = 0
                line_draw_step[0][4] = 0
                line_num += 1
                print("已选择元器件终点")
                return
            # 选中一般点
            elif check_point() == 0 or check_point() == 2:
                draw_line()
                line_draw_step[0][1] = line_draw_step[0][3]
                line_draw_step[0][2] = line_draw_step[0][4]

            elif check_point() == 4:
                print("error retry")
                line_draw_step[0][3] = 0
                line_draw_step[0][4] = 0
                return
        line_draw_step[0][0] = 1
        return

    if not selected:
        return
    # 获取鼠标点击的位置坐标
    # 计算鼠标点击位置对应的网格点的行号和列号（从0开始）
    row = (y - 50) // WIDTH
    col = (x - 50) // WIDTH
    # 如果行号或列号超出范围，则不做任何操作，直接返回
    if row < 0 or row >= 100 or col < 0 or col >= 100:
        return
    # 获取当前选中的元器件的类型和位置（如果有）
    name, old_row, old_col, image_id = selected
    # 如果鼠标点击位置与当前选中元器件的位置相同，则不做任何操作，直接返回（避免重复绘制）
    if row == old_row and col == old_col:
        return
    # 如果当前选中元器件已经在画布上绘制，则删除原来的图片
    if image_id is not None:
        canvas.delete(image_id)
    # 在画布上绘制新的元器件，并获取新的图片ID
    image_id = draw_image(name, row, col)

    # 更新当前选中元器件的位置和图片ID
    selected = [name, row, col, image_id]


def close_popup(window, entry1, entry2):
    # 获取输入框中的数值
    var.set(1)
    global value1, value2
    value1 = entry1.get()
    value2 = entry2.get()
    # 销毁窗口
    window.destroy()


def popup(name):
    window = tk.Toplevel()
    window.title("设置参数")
    window.geometry("230x150+400+200")
    window.config(bg="lightblue")

    if name == 'res':
        label1 = tk.Label(window, text="Ω：", bg="lightblue")
        label2 = tk.Label(window, text="？", bg="lightblue")
    if name == 'cap':
        label1 = tk.Label(window, text="F：", bg="lightblue")
        label2 = tk.Label(window, text="？", bg="lightblue")
    if name == 'ind':
        label1 = tk.Label(window, text="H：", bg="lightblue")
        label2 = tk.Label(window, text="？", bg="lightblue")
    if name == 'udc':
        label1 = tk.Label(window, text="U：", bg="lightblue")
        label2 = tk.Label(window, text="初相位", bg="lightblue")
    if name == 'idc':
        label1 = tk.Label(window, text="A：", bg="lightblue")
        label2 = tk.Label(window, text="初相位", bg="lightblue")
    entry1 = tk.Entry(window)
    entry2 = tk.Entry(window)

    button = tk.Button(window, text="确定",
                       command=lambda: close_popup(window, entry1, entry2))

    label1.grid(row=0, column=0, padx=10, pady=10)
    entry1.grid(row=0, column=1, padx=10, pady=10)
    label2.grid(row=1, column=0, padx=10, pady=10)
    entry2.grid(row=1, column=1, padx=10, pady=10)
    button.grid(row=2, column=1, padx=10, pady=10)
    button.wait_variable(var)

# close_popup_freq + popup_freq：频率的获取


def close_popup_freq(window, entry):
    var.set(1)
    global freq
    freq = entry.get()
    window.destroy()


def popup_freq():
    window = tk.Toplevel()
    window.title("设置频率")
    window.geometry("230x150+400+200")
    window.config(bg="lightblue")
    label = tk.Label(window, text="Hz", bg="lightblue")
    entry = tk.Entry(window)
    button = tk.Button(window, text="确定",
                       command=lambda: close_popup_freq(window, entry))
    label.grid(row=0, column=0, padx=10, pady=10)
    entry.grid(row=0, column=1, padx=10, pady=10)
    button.grid(row=1, column=1, padx=10, pady=10)
    button.wait_variable(var)

# 定义一个函数，用于处理键盘按键事件


def on_key(event):
    global selected  # 声明全局变量selected
    # 如果当前没有选中任何元器件，则不做任何操作，直接返回

    if not selected:
        return
    # 获取当前选中的元器件的类型和位置（如果有）
    name, row, col, image_id = selected
    # 根据按键的不同，更新元器件的位置
    if event.keysym == "Up":  # 则行号减一
        row -= 1
    elif event.keysym == "Down":  # 则行号加一
        row += 1
    elif event.keysym == "Left":  # 则列号减一
        col -= 1
    elif event.keysym == "Right":  # 则列号加一
        col += 1
    elif event.keysym == "space":  # 表示确认放置元器件，清空当前选中的元器件，并直接返回
        selected = []
        canvas.create_oval(50+(col-1)*WIDTH, 50+(row+3)*WIDTH, 50 +
                           (col-1)*WIDTH, 50+(row+3)*WIDTH, outline="red", width=5)
        canvas.create_oval(50+(col+1)*WIDTH, 50+(row+3)*WIDTH, 50 +
                           (col+1)*WIDTH, 50+(row+3)*WIDTH, outline="blue", width=5)

        popup(name)
        real_selected.append(
            (name, (col-1), (row+3), (col+1), (row+3), image_id, value1, value2))
        # print(real_selected)

        return
    else:  # 如果按下了其他键，则不做任何操作，直接返回
        return
    # 如果行号或列号超出范围，则不做任何操作，直接返回
    if row < 0 or row >= 100 or col < 0 or col >= 100:
        return
    # 计算元器件的新的中心坐标
    x = 50 + col * WIDTH
    y = 50 + row * WIDTH
    # 移动画布上的元器件图片到新的位置

    canvas.move(image_id, x - (50 + old_col * WIDTH),
                y - (50 + old_row * WIDTH))
    old_col = x
    old_row = y

    # 更新当前选中元器件的位置和图片ID
    selected = [name, row, col, image_id]


def draw_line():
    # print("in draw_line")
    canvas.create_line(50+line_draw_step[0][1]*WIDTH, 50+line_draw_step[0][2]
                       * WIDTH, 50+line_draw_step[0][3]*WIDTH, 50+line_draw_step[0][4]*WIDTH)


def get_matrix(matrix):
    x, y = np.shape(matrix)
    final_matrix = np.zeros((1, y), dtype=int)
    num = x
    while 1:
        node = np.zeros((1, y), dtype=int)
        for i in range(0, x):
            if not np.all(matrix[i, :].reshape(1, -1) == 0):
                node = matrix[i, :].reshape(1, -1).copy()
                matrix[i, :] = 0
                num -= 1
                break
        for i in range(0, x):
            tag = 0
            for j in range(0, y):
                if matrix[i][j] != 0 and matrix[i][j] == node[0][j]:
                    tag = 1
                    break
            if tag == 1:
                for j in range(0, y):
                    node[0][j] = node[0][j] if node[0][j] != 0 else matrix[i][j]
                matrix[i, :] = 0
                num -= 1
        final_matrix = np.concatenate((final_matrix, node), axis=0)
        if num == 0:
            break
    return final_matrix[1:]


def update_df(log):
    global incidence_matrix, df
    incidence_matrix = incidence_matrix[~(
        incidence_matrix == 0).all(axis=1)]  # 删除所有整行为0的行

    final_matrix = get_matrix(incidence_matrix.copy())
    # incidence_matrix：Version 1导线版
    # final_matrix: Version 2节点版

    if log:
        print("\033[1;33m________________Version 1________________\033[0m")
        print("连接表：（x:导线，y:元器件）\n可以理解为一根导线连接两个元器件，红色端电流流入记作（1），蓝色段电流流出记作（-1）\n",incidence_matrix)
        print("上述关联矩阵每一列对应一个元器件，分别是\n每一行各个元素为（元器件类型，x红，y红，x蓝，y蓝，器件ID）\n", real_selected)
        print("\033[1;32m________________Version 2________________\033[0m")
        print("关联矩阵：（x:节点，y：元器件）\n电流方向背离此节点记作（+1），流入此节点记作（-1），无关联记作0", final_matrix)
        print("上述关联矩阵每一列对应一个元器件，分别是\n（元器件类型，x红，y红，x蓝，y蓝，器件ID）\n", real_selected)

    df = pd.DataFrame(columns=["type", "start", "end", "value", "phi", "freq"])

    for j in range(len(final_matrix[0])):  # 元件
        this_list = ["" for i in range(6)]
        for i in range(len(final_matrix)):  # 节点
            this_list[0] = real_selected[j][0]
            if final_matrix[i][j] == -1:
                this_list[1] = i+1
            if final_matrix[i][j] == 1:
                this_list[2] = i+1
            this_list[3] = real_selected[j][6]
            this_list[4] = real_selected[j][7]
        df.loc[len(df)] = this_list

    df.loc[0, "freq"] = freq
    if log:
        print("\033[1;31m________________Version 3________________\033[0m")
        print("关联矩阵修改版：x:元器件；y:元器件类型(type)、起始节点(start)、终止结点(end)、值(value)、相位（电流源/电压源）(phi),频率记录在第0行的最后一个值(freq)")
        print(df)


# 定义一个函数，用于处理工具栏按钮点击事件


def on_button(name):
    global mode, incidence_matrix, df, selected
    # 绘制导线
    if name == "1":
        print("进入绘制导线模式，请勿重复点击按键")
        num = np.shape(real_selected)[0]
        incidence_matrix = np.zeros(shape=(num*num, num), dtype=int)
        mode = 1
        return

    # 设置频率(检验无误)
    if name == "2":
        popup_freq()
        print(freq)
        return

    # 查询关联矩阵 
    if name == "3":
        update_df(log=True)

    # 求解节点电压
    if name == "4":
        update_df(log=True)
        print("\033[1;31m_____________基于Version 3的求解模块_____________\033[0m")
        solve(df)

    # 如果当前已经选中了某个元器件，则删除原来的图片（如果有）
    if selected:
        _, _, _, image_id = selected
        if image_id is not None:
            canvas.delete(image_id)
    selected = [name, 0, 0, None]


def main():
    frame = tk.Frame(window)

    label = tk.Label(frame, text="工具栏", font=("Arial", 16))

    label.grid(row=0, column=0, sticky=tk.NSEW, padx=10, pady=10)
    # 存储工具栏按钮对应的元器件类型和名称
    buttons = [("电容", "cap"), ("电感", "ind"), ("电压源", "udc"), ("电流源", "idc"),
               ("电阻", "res"), ("绘制导线", "1"), ("设置频率", "2"), ("查询关联矩阵", "3"), ("求解", "4")]

    for i, (name, type) in enumerate(buttons):

        button = tk.Button(
            frame, text=name, command=lambda t=type: on_button(t))
        button.grid(row=i+1, column=0, sticky=tk.NSEW, padx=10, pady=10)

    frame.pack(side=tk.RIGHT)
    canvas.pack(side=tk.LEFT)
    # 绑定鼠标左键点击事件到画布对象上，调用on_click函数
    canvas.bind("<Button-1>", on_click)
    # 绑定键盘按键事件到窗口对象上，调用on_key函数
    window.bind("<Key>", on_key)
    # 调用draw_grid函数，在画布上绘制网格点
    draw_grid()
    # 启动窗口对象的主循环
    window.mainloop()


if __name__ == "__main__":
    main()
