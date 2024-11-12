import math

def error(txt):
    raise ValueError("\033[31m" + txt + "\033[39m") 

def evaluate(value, reference, threshold, operator = "relative",  txt = ""):

    flag = False
    error_value = 0

    if (operator == "<"):

        error_value = math.fabs(value - threshold)

        flag = (value > threshold)

        if (flag):
            error(txt + ": {} > {} with error {}".format(value, threshold, error_value)) 

    elif (operator == "=="):

        error_value = math.fabs(value - threshold)

        flag = (error_value != 0)

        if (flag):
            error(txt + ": {} not equal to {} with error {}".format(value, threshold, error_value)) 

    elif (operator == "relative"):

        if (reference == 0):
            error("Can not evaluate a relative error with reference == 0")

        error_value = math.fabs((value - reference)/reference)

        flag = (error_value > threshold)

        if (flag):
            error(txt + ": {} vs {} with error {} for threshold {}".format(value, reference, error_value, threshold)) 

    elif (operator == "absolute"):

        error_value = math.fabs(value - reference)

        flag = (error_value > threshold)

        if (flag):
            error(txt + ": {} vs {} with error {} for threshold {}".format(value, reference, error_value, threshold)) 

    else:

        error("Operator not recognized")
