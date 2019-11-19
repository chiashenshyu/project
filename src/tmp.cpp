#include "common.hpp"
#include "ParticleFilter.hpp"
#include "planner.hpp"
#include "mpc_steering.hpp"

#define VIZ

int calc_lookahead_pt(std::vector<double> cx,std::vector<double> cy,States veh){
    double min_dis = std::numeric_limits<double>::max();
    int indx = 0, i = 0;
    for(;i<cx.size();i++){
        double dx = veh.x - cx[i];
        double dy = veh.y - cy[i];
        double tempdist = sqrt(pow(dx,2)+pow(dy,2));
        if(min_dis>tempdist){
            min_dis = tempdist;
            indx = i;
        }
    }
    // double angle = pi2pi(cd[indx] - atan2(cy[indx] - veh.y, cx[indx] - veh.x));
    // // if(angle < 0) indx = -indx;
    // isRight = (angle < 0)? true : false; 
    return indx;
}

void get_14points(int indx,Eigen::VectorXd &x_14pts,Eigen::VectorXd &y_14pts,
                  std::vector<double>cx,std::vector<double>cy,std::vector<double>cd,
                  double lookaheadDist){
    // for(int i=0;i<fit_numberof_pts;i++){
    //     x_14pts[i] = cx[i+indx];
    //     y_14pts[i] = cy[i+indx];
    // }
    int i = 0;
    double totalDist = 0.0;  
    std::vector<double> x, y; 
    while((totalDist < lookaheadDist || i < 14) && i+indx < cx.size()){
        totalDist += (i != 0)? cd[i+indx] : 0;
        i++;
    }
    x_14pts.resize(i);
    y_14pts.resize(i); 
    for(int j = 0; j < i; j++){
        x_14pts[j] = cx[j+indx];
        y_14pts[j] = cy[j+indx];
    }
    // cout << "points selected: " << i << endl;
}

void transform_to_local(States veh_,Eigen::VectorXd &x_14pts,Eigen::VectorXd &y_14pts){
    Eigen::VectorXd pnt(2);
    Eigen::VectorXd local_pnt(2);

    Eigen::MatrixXd translation(2,2);
    translation <<  cos(-veh_.theta), -sin(-veh_.theta),
                    sin(-veh_.theta),  cos(-veh_.theta);
    for(int i =0; i<x_14pts.size();i++){
        pnt << x_14pts[i] - veh_.x, y_14pts[i] - veh_.y;
        local_pnt = translation * pnt;
        x_14pts[i] = local_pnt[0];
        y_14pts[i] = local_pnt[1];
    }
}

void transform_to_global(States veh, vector<double>& x, vector<double>& y){
    assert(x.size() == y.size()); 
    int n = x.size();
    double theta = veh.theta, xtmp, ytmp; 
    for(int i = 0; i < n; i++){
        xtmp = x[i]; 
        ytmp = y[i];
        x[i] = veh.x + xtmp*cos(theta) - ytmp*sin(theta); 
        y[i] = veh.y + xtmp*sin(theta) + ytmp*cos(theta); 
    }
}

void transform_to_local(States veh, std::vector<std::vector<double>>& obstacles,
                        std::vector<std::vector<double>>& obstacles_local){
    for(int i = 0; i < obstacles.size(); i++){
        double xtmp = obstacles[i][0] - veh.x, ytmp = obstacles[i][1] - veh.y; 
        obstacles_local[i][0] = cos(-veh.theta) * xtmp - sin(-veh.theta) * ytmp; 
        obstacles_local[i][1] = sin(-veh.theta) * xtmp + cos(-veh.theta) * ytmp;
    }
}

double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

vector<double> polyvec (Eigen::VectorXd coeffs, vector<double>& x,int points){
    vector<double> y;
    double x_ = x[0];
    y.push_back(polyeval(coeffs,x_));
    for(int i = 1;i<points;i++){
        y.push_back(polyeval(coeffs,x_));
        x.push_back(x_);
        x_ += 0.1;
    }
    return y;
}

void polyfit(Eigen::VectorXd& coeff, Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
    assert(xvals.size() == yvals.size());
    int err = -1; 
    if(!(order >= 1 && order <= xvals.size() - 1)){
        throw err; 
    }
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    coeff = Q.solve(yvals);
    return;
}

int main(){
    /**
     * RRT star planner
     */
    planner_params A;
    A.origin = Point( 400, -400); // start point
    A.goal   = Point(-400,  400); // end point

    // Maze Map
    MatrixXd obstacleMap(5,4); 
    obstacleMap.row(0) << 50,0,50,150;
    obstacleMap.row(1) << 150,200,150,100;
    obstacleMap.row(2) << 150,100,100,100;
    obstacleMap.row(3) << 125,0,125,50;
    obstacleMap.row(4) << 200,75,150,75;
    A.obstacle = obstacleMap - (100*MatrixXd::Ones(5,4));
    A.obstacle.col(0) = -1*A.obstacle.col(0);
    A.obstacle.col(2) = -1*A.obstacle.col(2);
    A.obstacle *= 5; 

    // Map params
    A.iterations = 40000;
    A.width      = 1000; 
    A.height     = 1000;
    A.goalProx   = 15;

    // RRT planner
    Planner p(A); 
    if(p.RRTstar()){
        std::cout << "Failed to find path to goal point" << std::endl;
        return 1; 
    }

    // Ectract path from planner for controller 
    Path path; 
    std::vector<Node> wayPoints;
    std::vector<double> _x, _y; 
    p.ExtractPath(path, wayPoints);
    std::reverse(path.cx.begin(), path.cy.end()); 
    std::reverse(path.cy.begin(), path.cy.end()); 
    smooth(path.cx, path.cy, _x, _y); 
    path.cx = _x; 
    path.cy = _y;
    std::vector<double> cd(path.cx.size(), 0); 
    for(int i = 1; path.cy.size(); i++){
        cd[i] = sqrt(0.01 + pow(path.cy[i]-path.cy[i-1], 2)); 
    }

    /**
     * Visualizer
     */
    #ifdef VIZ 
    Visualizer viz; 
    viz.plannerParamsIn(A);
    #endif

    /**
     * PF Localization
     */
    ParticleFilter pf(1000); 
    std::vector<std::vector<double>> landmark = {{   0,  400},
                                                 {-450,  450},
                                                 {-400, -400},
                                                 { 100, -300},
                                                 { 500,  350},
                                                 { 425, -450},
                                                 {-155,  200},
                                                 {-550,   50},
                                                 { 200,  530},
                                                 { 550,    0},
                                                 { -50, -475}
                                                 };
    std::vector<bool> lmVisible(landmark.size(), false);
    std::vector<std::vector<double>> errorArray(landmark.size()+1);

    /**
     *  MPC controller
     */
    MPC mpcController; 

    // Obstacles
    std::vector<std::vector<double>> obstacles, obstaclesLocal; 
    // obstacles => todo 
    obstaclesLocal = obstacles; 

    // Set initial condition
    States veh(0., 3., 0., 2.);
    std::vector<States> vehStateVec; 
    std::vector<double> x, y, v, t; 
    vehStateVec.push_back(veh); 

    // Params 
    double maxTime = 1000., _time = 0., delta = 0., a = 1., lookaheadDist = 0.; 
    int fittingPoints = 14, lastTargetIndex = 0; 
    Eigen::VectorXd x_14pts(fittingPoints); 
    Eigen::VectorXd y_14pts(fittingPoints); 

    while(_time < maxTime){
        // Take 14 points from the path, starting from the closest point and transform 
        // them into local frame
        transform_to_local(veh, obstacles, obstaclesLocal);
        int closetIndex = calc_lookahead_pt(path.cx, path.cy, veh);
        get_14points(closetIndex, x_14pts, y_14pts, path.cx, path.cy, cd, lookaheadDist); 
        std::vector<double> x14pts(x_14pts.size(), 0.), y14pts(y_14pts.size(), 0.); 
        for(int i = 0; i < x_14pts.size(); i++){
            x14pts[i] = x_14pts[i]; 
            y14pts[i] = y_14pts[i]; // check this part, seems to be unnecessary
        }
        transform_to_local(veh, x_14pts, y_14pts); 
        
        // Fit a 3 order spline to those chosen points
        Eigen::VectorXd coeff; 
        try{
            polyfit(coeff, x_14pts, y_14pts, 3); 
        }catch(int e){
            std::cout << "Done!" << std::endl;
            break; 
        }

        // Set initial state for MPC controller
        double cte = coeff[0], we = -atan(coeff[1]); 
        Eigen::VectorXd vehState(6); 
        vehState << 0., 0., 0., veh.v, cte, we; 

        // Run MPC controller
        auto res = mpcController.mpc_solve(vehState, coeff, obstaclesLocal); 
        veh.update(res[1], res[0], 0.1);

        // Set lookahead distance for next step according to vehicle state now
        double xtmp = 0., ytmp = 0.;
        std::vector<double> xmpc, ympc; 
        for(int i = 2; i < res.size(); i+= 2){
            lookaheadDist += sqrt(pow(res[i]-xtmp,2) + pow(res[i+1]-ytmp,2)); 
            xtmp = res[i]; 
            ytmp = res[i+1]; 
            xmpc.push_back(xtmp); 
            ympc.push_back(ytmp);
        }
        transform_to_global(veh, xmpc, ympc); 

        _time += 0.1; 
        vehStateVec.push_back(veh);
        lastTargetIndex = closetIndex; // Fix this in the future
        
        /**
         * Visualization
         */
        #ifdef VIZ
            x.push_back(vehStateVec.back().x); 
            y.push_back(vehStateVec.back().y);
            // Clear previous plot
            plt::clf(); 
            // Plot line from given x and y data. Color is selected automatically. 
            plt::named_plot("Car", x, y, "*k");
            plt::named_plot("Path", path.cx, path.cy, "-r");
            plt::named_plot("mpc_traj", xmpc, ympc, "go");
            for(int i = 0; i < obstacles.size(); i++) plotCircle(obstacles[i][0], obstacles[i][1], 1.0);
            // Legends on; Axis equal
            plt::legend(); 
            plt::axis("equal");
            plt::pause(0.1);
        #endif
        }
        plt::show();     
    return 0;
}