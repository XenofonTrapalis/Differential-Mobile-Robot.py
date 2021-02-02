#Python 3.6.4 (v3.6.4:d48eceb, Dec 19 2017, 06:54:40) [MSC v.1900 64 bit (AMD64)] on win32
#Type "copyright", "credits" or "license()" for more information.
# https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html
import numpy as np
import pandas as pd
from numpy.linalg import inv
from math import *
from array import *
import matplotlib.pyplot as plt
import os

#import data from csv file
data = pd.read_csv(r'..\Your Path\import_data\differential mobile robot data.csv', delimiter = ',' )
#import data from User
radius = float(input("Δώσε ακτίνα τροχών se cm" + '\n'))/100 
baseline = float(input("Δώσε την απόσταση από την μέση του άξονα του τροχού se cm" + '\n'))/100
#lists for data manipulation
data_list = [] 
data_time_sec = [] 
data_rpm_left_df_dt = [] 
data_rpm_right_df_dt = []
data_rpm_left_df_dt_new = []
data_rpm_right_df_dt_new = []
data_time_sec_new = []
#constants for start point
theta = 0
vx = 0
vy = 0
vtheta = 0
#wheel constants
ar = -pi/2         #right wheel
br = 0
lr = baseline
al = -pi/2         #left wheel
bl = 0
ll = -baseline
temp = np.zeros((3,1)) 

usvcc = np.array([[vx], [vy], [vtheta]]) #pose array

Rolling_constrains = np.array([[-sin(ar+br), cos(ar+br), lr*cos(br)],  #rolling constrains array
                               [-sin(al+bl), cos(al+bl), ll*cos(bl)]])

No_sliding_constrains = np.array([[0, -1, 0]]) #no sliding constrain array

radius_array = np.array([[radius, 0],   #radius array
                         [0, radius]])


#manipulation data from csv file
for key in data:
    data_list.append(data[key])
for element in data_list[0]:
    data_time_sec.append(element)
for element in data_list[1]:
    data_rpm_left_df_dt.append(element)
for element in data_list[2]:
    data_rpm_right_df_dt.append(element)
for index,element in enumerate(data_rpm_left_df_dt):
    data_rpm_left_df_dt[index] = element * 2*pi/60
for index,element in enumerate(data_rpm_right_df_dt):
    data_rpm_right_df_dt[index] = element * 2*pi/60
#fill in the missing times and data
index = 30
for time_p,time_n,sppedL,speedR in zip(data_time_sec,data_time_sec[1:],data_rpm_left_df_dt,data_rpm_right_df_dt):
    if abs(time_p - time_n)> 1:
        w = abs(time_p - time_n)
        for w in range(index,index+w):
            data_time_sec_new.insert(index,time_p)
            data_rpm_left_df_dt_new.insert(index,sppedL)
            data_rpm_right_df_dt_new.insert(index,speedR)
    else:
        index += 1
        data_time_sec_new.insert(index,time_p)
        data_rpm_left_df_dt_new.insert(index,sppedL)
        data_rpm_right_df_dt_new.insert(index,speedR)       
data_time_sec_new.append(data_time_sec[-1])
data_rpm_left_df_dt_new.append(data_rpm_left_df_dt[-1])
data_rpm_right_df_dt_new.append(data_rpm_right_df_dt[-1])
for index,time in enumerate(data_time_sec_new):
    data_time_sec_new[index] = index
print("data" + '\n' + "sec" + "|" + "right wheel rad/sec" + "|" + "left wheel rad/sec")
print(*zip(data_time_sec_new,data_rpm_right_df_dt_new,data_rpm_left_df_dt_new), sep='\n')
#constrain and velocity array
speed_of_wheel_array = np.vstack((data_rpm_right_df_dt_new,data_rpm_left_df_dt_new))
speed_of_wheel_array = np.matmul(radius_array,speed_of_wheel_array)
speed_of_wheel_array = np.vstack((speed_of_wheel_array,np.zeros((1,31))))
constrains_array = np.vstack((Rolling_constrains,No_sliding_constrains))
#lists for export data to csv file
robot_data_x= []
robot_data_y = []
robot_data_theta = []
global_references_data_x = []
global_references_data_y = []
global_references_data_theta = []
pose_data_x = []
pose_data_y = []
pose_data_theta = []
taxythes_adraneiako_systhma_data_right = []
taxythes_adraneiako_systhma_data_left = []
taxythtes_systhma_robot_data_right = []
taxythtes_systhma_robot_data_left = []
rpm_data_right = []
rpm_data_left = []
#kinhmatiko kai antistrofo kinhmatiko provlhma
for index in range(speed_of_wheel_array.shape[1]):
#kinhmatiko provlhma
    robot_references = np.matmul(inv(constrains_array),speed_of_wheel_array[:3 , index:index+1])
    Rtheta = np.array([[cos(theta), sin(theta), 0],
                       [-sin(theta), cos(theta), 0],
                       [0, 0, 1]])
    global_references = np.matmul(inv(Rtheta),robot_references)
    theta = global_references[2] + theta
    usvcc = np.add(global_references,usvcc)
    plt.plot(usvcc[0],usvcc[1], 'bo')
#export data manipulation
    robot_data_x.append(robot_references[0].tolist())
    robot_data_y.append(robot_references[1].tolist())
    robot_data_theta.append(robot_references[2].tolist())
    global_references_data_x.append(global_references[0].tolist())
    global_references_data_y.append(global_references[1].tolist())
    global_references_data_theta.append(global_references[2].tolist())
    pose_data_x.append(usvcc[0].tolist())
    pose_data_y.append(usvcc[1].tolist())
    pose_data_theta.append(usvcc[2].tolist())
#end of export data manipulation
#antistrofo kinhmatiko provlhma
    global_references_anti_kin = np.subtract(usvcc,temp)
    temp = usvcc
    Rtheta_ant_kin = np.array([[cos(temp[2]), sin(temp[2]), 0],
                       [-sin(temp[2]), cos(temp[2]), 0],
                       [0, 0, 1]])
    a = np.matmul(inv(radius_array),Rolling_constrains)
    b = np.matmul(Rtheta,global_references_anti_kin)
    taxythtes_systhma_robot = np.matmul(a,global_references_anti_kin)
    taxythes_adraneiako_systhma = np.matmul(a,b)
    rpm = taxythes_adraneiako_systhma * 60./(2.*pi)
#export data manipulation
    taxythes_adraneiako_systhma_data_right.append(taxythes_adraneiako_systhma[0].tolist())
    taxythes_adraneiako_systhma_data_left.append(taxythes_adraneiako_systhma[1].tolist())
    taxythtes_systhma_robot_data_right.append(taxythtes_systhma_robot[0].tolist())
    taxythtes_systhma_robot_data_left.append(taxythtes_systhma_robot[1].tolist())
    rpm_data_right.append(rpm[0].tolist())
    rpm_data_left.append(rpm[1].tolist())
#end of export data manipulation
    


if not os.path.isdir('./Export_data'):
    directory = "Export_data"
    parent_dir = "..\Trapalis_project_1"
    path = os.path.join(parent_dir,directory)
    os.mkdir(path)
robot_references_data = {'sec':[data_time_sec_new],'x':[robot_data_x],'y':[robot_data_y],'theta':[robot_data_theta]}
robot_data = pd.DataFrame(robot_references_data, columns=['sec','x','y','theta'])
robot_data.to_csv (r'..\Your Path\Export_data\dianisma_taxythtwn_tou_robot_data.csv', index = False, header=True)

global_references_data = {'sec':[data_time_sec_new],'x':[global_references_data_x],'y':[global_references_data_y],'theta':[global_references_data_theta]}
global_data = pd.DataFrame(global_references_data, columns=['sec','x','y','theta'])
global_data.to_csv (r'..\Your Path\Export_data\dianisma_taxythtwn_adraneiko_systhma_data.csv', index = False, header=True)

pose_data = {'sec':[data_time_sec_new],'x':[pose_data_x],'y':[pose_data_y],'theta':[pose_data_theta]}
pose = pd.DataFrame(pose_data, columns=['sec','x','y','theta'])
pose.to_csv (r'..\Your Path\Export_data\pose_data.csv', index = False, header=True)

taxythes_adraneiako_systhma_data = {'sec':[data_time_sec_new],'right wheel rad/sec':[taxythes_adraneiako_systhma_data_right],'left wheel rad/sec':[taxythes_adraneiako_systhma_data_left]}
taxythes_adraneiako_data = pd.DataFrame(taxythes_adraneiako_systhma_data, columns=['sec','right wheel rad/sec','left wheel rad/sec'])
taxythes_adraneiako_data.to_csv (r'..\Your Path\Export_data\taxythes_adraneiako_systhma_data.csv', index = False, header=True)

taxythtes_robot_data = {'sec':[data_time_sec_new],'right wheel rad/sec':[taxythtes_systhma_robot_data_right],'left wheel rad/sec':[taxythtes_systhma_robot_data_left]}
taxythtes_systhma_robot_data = pd.DataFrame(taxythtes_robot_data, columns=['sec','right wheel rad/sec','left wheel rad/sec'])
taxythtes_systhma_robot_data.to_csv (r'..\Your Path\Export_data\taxythtes_systhma_robot_data.csv', index = False, header=True)

rpm_data_wheel = {'sec':[data_time_sec_new],'right wheel df/dt':[rpm_data_right],'left wheel df/dt':[rpm_data_left]}
rpm_data = pd.DataFrame(rpm_data_wheel, columns=['sec','right wheel df/dt','left wheel df/dt'])
rpm_data.to_csv (r'..\Your Path\Export_data\rpm_data.csv', index = False, header=True)

print('\n' + '\n' + "kinhmatiko provlhma",'\n')
print("dianisma taxuthtwn tou robot os pros to systhma R tou robot",'\n')
print("sec" + " | " + " x rad/sec" + " | " + " y rad/sec" + " | " +" theta rad/sec")
print(*zip(data_time_sec_new,robot_data_x,robot_data_y,robot_data_theta),sep='\n')
print('\n' + '\n' + "dianisma taxuthtwn os pros to adraneiako systhma",'\n')
print("sec" + " | " + " x rad/sec" + " | " + " y rad/sec" + " | " + " theta rad/sec")
print(*zip(data_time_sec_new,global_references_data_x,global_references_data_y,global_references_data_theta),sep='\n')
print('\n' + '\n' + "thesh tou robot sto systhma tou edafous",'\n')
print("sec" + " | " + " x rad/sec" + " | " + " y rad/sec" + " | " + " theta rad/sec")
print(*zip(data_time_sec_new,pose_data_x,pose_data_y,pose_data_theta),sep='\n')

print('\n' + '\n' + "antistrofo kinhmatiko provlhma",'\n')
print("taxythtes os pros to adraneiako systhma",'\n')
print("sec" + " | " + " df/dt right in rpm" + " | " + " df/dt left in rpm")
print(*zip(data_time_sec_new,taxythes_adraneiako_systhma_data_right,taxythes_adraneiako_systhma_data_left),sep='\n')
print('\n' + '\n' + "taxythtes sto systhma anaforas R tou robot",'\n')
print("sec" + " | " + " df/dt right in rpm" + " | " + " df/dt left in rpm")
print(*zip(data_time_sec_new,taxythtes_systhma_robot_data_right,taxythtes_systhma_robot_data_left),sep='\n')
print('\n' + '\n' + "dianysma taxythtwn df/dt se rpm",'\n')
print("sec" + " | " + " df/dt right in rpm" + " | " + " df/dt left in rpm")
print(*zip(data_time_sec_new,rpm_data_right,rpm_data_left),sep='\n')

plt.show()
