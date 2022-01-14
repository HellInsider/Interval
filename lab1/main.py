from functions import Function, ThreeHumpCamelFunction, HimmelblauFunction, SphereFunction, BoothFunction
from globopt0 import Globopt


def run(function: Function, eps: float) -> None:
    globopt = Globopt(function)
    est, work_list = globopt.globopt0(eps)

    print('Lead boxes')
    for i in range(5):
        print(work_list.at(i).domain.center())

    globopt.plot_boxes(work_list)
    globopt.plot_boxes_centers(work_list)
    globopt.plot_boxes_rads(work_list)
    globopt.plot_convergence(work_list)


def main():
    #run(HimmelblauFunction(), 0.00001)
    run(BoothFunction(), 0.00001)
    #run(ThreeHumpCamelFunction(), 0.00001)
    return


if __name__ == '__main__':
    main()
