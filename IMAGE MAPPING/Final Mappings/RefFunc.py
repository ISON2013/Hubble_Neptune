import numpy as np
import matplotlib.pyplot as plt

def res_to_lat(res,lat,val):
    ratio = (np.max(lat)-np.min(lat))/res
    val = val*ratio
    return val

def lat_to_res(res,lat,val):
    ratio = res/(np.max(lat)-np.min(lat))
    val = val*ratio
    return val

def res_to_long(res,long,val):
    ratio = (np.max(long)-np.min(long))/res
    val = val*ratio
    return val

def long_to_res(res,long,val):
    ratio = res/(np.max(long)-np.min(long))
    val = val*ratio
    return val

def restrict_long(img,right,left):
    if right == '':
        right = img.shape[1]
    if left == '':
        left = 0
    right = int(right)
    left = int(left)
    img = img[0:img.shape[0],left:right]
    return img