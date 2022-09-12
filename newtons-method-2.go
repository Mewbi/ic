package main

import (
	"fmt"
	"math"
)

var A [4]float64 = [4]float64{-200, -100, -170, 15}
var a [4]float64 = [4]float64{-1, -1, -6.5, 0.7}
var b [4]float64 = [4]float64{0, 0, 11, 0.6}
var c [4]float64 = [4]float64{-10, -10, -6.5, 0.7}
var x_0 [4]float64 = [4]float64{1, 0, -0.5, -1}
var y_0 [4]float64 = [4]float64{0, 0.5, 1.5, 1}

type grad func(float64, float64) [2]float64
type hess func(float64, float64) [4]float64

func inv_matrix(A [4]float64) [4]float64 {
    det := float64(A[0]*A[3] - A[1]*A[2])
    return [4]float64{A[3]/det, -A[1]/det, -A[2]/det, A[0]/det}
}


func mult_matrix(A [4]float64, v [2]float64) [2]float64{
    return [2]float64{A[0]*v[0] + A[2]*v[1], A[1]*v[0] + A[3]*v[1]}
}

func f(i int, x, y float64) float64 {
    aux_x := x - x_0[i]
    aux_y := y - y_0[i]
    return a[i]*aux_x*aux_x + b[i]*aux_x*aux_y + c[i]*aux_y*aux_y
}

func partial_f_x(i int, x, y float64) float64 {
    aux_x := x - x_0[i]
    aux_y := y - y_0[i]
    return 2*a[i]*aux_x + b[i]*aux_y
}

func partial_f_y(i int, x, y float64) float64 {
    aux_x := x - x_0[i]
    aux_y := y - y_0[i]
    return b[i]*aux_x + 2*c[i]*aux_y
}

func grad_V(x, y float64) [2]float64 {
    var D_x, D_y float64
    for i := 0; i < 4; i++  {
        D_x += A[i]*math.Exp(f(i, x, y))*partial_f_x(i, x, y)
        D_y += A[i]*math.Exp(f(i, x, y))*partial_f_y(i, x, y)
    }
    return [2]float64{D_x , D_y}
}

func hess_V(x, y float64) [4]float64 {
    var D_xx, D_xy, D_yy float64
    for i := 0; i < 4; i++  {
        aux_x := partial_f_x(i, x, y)
        aux_y := partial_f_y(i, x, y)
        D_xx += A[i]*math.Exp(f(i, x, y))*(aux_x*aux_x + 2*a[i])
        D_yy += A[i]*math.Exp(f(i, x, y))*(aux_y*aux_y + 2*c[i])
        D_xy += A[i]*math.Exp(f(i, x, y))*(aux_x*aux_y + b[i])
    }
    return [4]float64{D_xx, D_xy, D_xy, D_yy}
}


func newton_step(f grad, fp hess, x, y float64) [2]float64{
    step := mult_matrix(inv_matrix(fp(x ,y)), f(x, y))
    return [2]float64{x - step[0], y - step[1]}
}


func newtons_method(f grad, fp hess, start [2]float64, n_iter int, tolerance float64) ([2]float64, error){
    var converge [2]float64
    prev_x := start[0]
    prev_y := start[1]
    for i := 0; i < n_iter; i++  {
        new_values := newton_step(f, fp, prev_x, prev_y)
        new_x := new_values[0]
        new_y := new_values[1]
        grad_values := f(new_x ,new_y)
        grad_x := grad_values[0]
        grad_y := grad_values[1]

        if math.Abs(grad_x) + math.Abs(grad_y) < tolerance {
            return [2]float64{new_x, new_y}, nil
        }
        prev_x = new_x
        prev_y = new_y
    }
    return converge, fmt.Errorf("Newton's Method didn't converge.")
}

func main() {
    c, err := newtons_method(grad_V, hess_V, [2]float64{0, 0}, 10000000, 0.0000000001)
    if err != nil {
        fmt.Println(err)
        return
    }
    fmt.Printf("Converge in %v \n", c)
}
