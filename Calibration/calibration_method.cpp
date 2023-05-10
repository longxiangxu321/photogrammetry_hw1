/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,    /// output: focal length (i.e., K[0][0], which is equal to 'alpha' in our slides).
        double& fy,    /// output: focal length (i.e., K[1][1], which is equal to 'beta/sin(theta)' in our slides).
        double& cx,    /// output: x component of the principal point (i.e., K[0][2], which is 'u0' in our slides).
        double& cy,    /// output: y component of the principal point (i.e., K[1][2], which is 'v0' in our slides).
        double& skew,  /// output: skew factor (i.e., K[0][1], which is equal to '-alpha * cot(theta)' in our slides).
        Matrix33& R,   /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t)   /// output：a 3D vector encoding camera translation.
{


//    Construct P

    int num_of_pairs = points_3d.size();
    Matrix P(2*num_of_pairs, 12);

    for (int i = 0; i < points_3d.size(); i++) {
        auto p3d = points_3d[i].homogeneous();
        auto up = - points_2d[i].x() * p3d;
        auto vp = - points_2d[i].y() * p3d;

        P.set_row(2*i, {p3d.x(), p3d.y(), p3d.z(), 1, 0, 0, 0, 0, up[0], up[1], up[2], up[3]});
        P.set_row(2*i+1, {0,0,0,0, p3d.x(), p3d.y(), p3d.z(), 1, vp[0], vp[1], vp[2], vp[3]});
    }

    double x=0.0;
    Matrix U(2*num_of_pairs, 2*num_of_pairs, x);
    Matrix D(2*num_of_pairs, 12, x);
    Matrix V(12,12,x);
    svd_decompose(P, U, D, V);

    Vector m = V.get_column(11);
    Vector3D a1({m[0], m[1], m[2]});
    Vector3D a2({m[4], m[5], m[6]});
    Vector3D a3({m[8], m[9], m[10]});
    double b1 = m[3];
    double b2 = m[7];
    double b3 = m[11];


    double rho, u0, v0, cos_theta, alpha, beta;

    rho = 1/a3.length();
    u0 = rho*rho*(dot(a1,a3));
    v0 = rho*rho*(dot(a2,a3));
    cos_theta = dot(cross(a1, a3), cross(a2,a3)) / ((cross(a1,a3).length() * cross(a2, a3).length()));
    double sin_theta = sqrt(1- cos_theta*cos_theta);
    alpha = rho*rho* cross(a1,a3).length()* sin_theta;
    beta = rho*rho* cross(a2, a3).length() * sin_theta;

    auto r1 = cross(a2, a3) / cross(a2,a3).length();
    auto r3 = rho * a3;
    auto r2 = cross(r3, r1) * 1;


//    double& fx,    /// output: focal length (i.e., K[0][0], which is equal to 'alpha' in our slides).
//    double& fy,    /// output: focal length (i.e., K[1][1], which is equal to 'beta/sin(theta)' in our slides).
//    double& cx,    /// output: x component of the principal point (i.e., K[0][2], which is 'u0' in our slides).
//    double& cy,    /// output: y component of the principal point (i.e., K[1][2], which is 'v0' in our slides).
//    double& skew,  /// output: skew factor (i.e., K[0][1], which is equal to '-alpha * cot(theta)' in our slides).
//    Matrix33& R,   /// output: the 3x3 rotation matrix encoding camera rotation.
//    Vector3D& t)   /// output：a 3D vector encoding camera translation.
    fx = alpha;
    fy = beta/sin_theta;
    cx = u0;
    cy = v0;
    skew = (-alpha) * (cos_theta/sin_theta);
    R.set_row(0,r1);
    R.set_row(1,r2);
    R.set_row(2,r3);

    Matrix33 K;
    K.set_row(0, {alpha, skew, cx});
    K.set_row(1, {0, fy, cy});
    K.set_row(2, {0,0,1});

    Vector3D b(b1,b2,b3);
    t = rho* inverse(K) * b;

    Matrix34 transform;
    Vector4D vv1(r1[0], r1[1], r1[2], t.x());
    Vector4D vv2(r2[0], r2[1], r2[2], t.y());
    Vector4D vv3(r3[0], r3[1], r3[2], t.z());
    transform.set_row(0, vv1);
    transform.set_row(1, vv2);
    transform.set_row(2, vv3);

    auto result1 = mult(transform, points_3d[5].homogeneous());
    auto result2 = mult(K, result1);
    auto result3 = result2/result2[2];
    std::cout<<"rho "<<rho<<std::endl;
    std::cout<<"input 3D "<<points_3d[5]<<std::endl;
    std::cout<<"result "<<result3<<std::endl;



    return true;
    return false;
}

















