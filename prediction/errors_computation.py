import numpy as np



def alt_error(alt_p, alt_t):
    return alt_p - alt_t


def ate_cte_error(x_p, y_p, trk_p, x_t, y_t):
    """ p - is prediction
        t - is true flight (radar track data)"""
    # this subtraction gives - for ate ahead and - for cte on the right
    d_x = x_p - x_t
    d_y = y_p - y_t

    # this gives + for ate ahead and + for cte on the right
    # d_x = x_p - x_b
    # d_y = y_p - y_b

    rad_trk = np.radians(trk_p)

    ate = d_x * np.sin(rad_trk) + d_y * np.cos(rad_trk)
    cte = d_x * np.cos(rad_trk) - d_y * np.sin(rad_trk)

    return ate, cte



if __name__ == "__main__":
    from math import cos, sin, radians

    # test 1
    x_b, y_b = np.array([10, 10, 10, 10, 10, -1]), np.array([10, 10, 10, 10, 10, -25])
    x_p, y_p = np.array([10, 10, 10, 10, 10, 30]), np.array([20, 20, 20, 20, 20, -40])

    trk_b = np.array([45, 0, 90, 225, 359, 60])
    print('\n', trk_b)
    print(ate_cte_error(x_b, y_b, trk_b, x_p, y_p))

    # trk_b = 0
    # print('\n', trk_b)
    # print(ate_cte_error(x_b, y_b, trk_b, x_p, y_p))
    #
    # trk_b = 90
    # print('\n', trk_b)
    # print(ate_cte_error(x_b, y_b, trk_b, x_p, y_p))
    #
    # trk_b = 225
    # print('\n', trk_b)
    # print(ate_cte_error(x_b, y_b, trk_b, x_p, y_p))
    #
    # trk_b = 359
    # print('\n', trk_b)
    # print(ate_cte_error(x_b, y_b, trk_b, x_p, y_p))