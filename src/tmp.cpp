#include "common.hpp"
#include "ParticleFilter.hpp"
#include "planner.hpp"
#include "mpc_steering.hpp"

#define VIZ

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
    Planner p(a); 
    if(p.RRTstar()){
        std::cout << "Failed to find path to goal point" << std::endl;
        return 1; 
    }

    // Ectract path from planner for controller 
    Path path; 
    std::vector<Node> wayPonts, _x, _y; 
    p.ExtractPath(path, wayPoints);
    std::reverse(path.cx.begin(), path.cy.end()); 
    std::reverse(path.cy.begin(), path.cy.end()); 
    smooth(path.cx, path.cy, _x, _y); 
    path.cx = _x; 
    path.cy = _y;
    std::vector<double> cd(path.cx.size(), 0); 
    for(int i = 1; cy.size(); i++){
        cd[i] = sqrt(0.01 + pow(cy[i]-cy[i-1], 2)); 
    }

    /**
     * Visualizer
     */
    #ifdef VIZ 
    Visualizer viz; 
    viz.plannerParamsIn(A): 
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
            y14pts[i] = y_14ptd[i]; // check this part, seems to be unnecessary
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
        auto res = controller.mpc_solve(vehState, coeff, obstaclesLocal); 
        veh.update(res[1], res[0], dt);

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

        time += dt; 
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
        }
        plt::show(); 
        #endif
    }
    
    return 0;
}