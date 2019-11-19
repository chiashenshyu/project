#include"mpc_steering.hpp"

using CppAD::AD;

int N = 20;
size_t x_start = 0;
size_t y_start = x_start+N;
size_t w_start = y_start+N;
size_t v_start = w_start+N;
size_t cte_start = v_start+N;
size_t we_start = cte_start+N;
size_t obs_start = we_start+N;
size_t delta_start = we_start+N;
size_t a_start = delta_start+N-1;
double target_speed = 5.0; //m/s
const double Lf = 2.67;
double dt = 0.1;
double _dx = 0.1; 



class FG_eval{
    public:
    Eigen::VectorXd coeff;
    std::vector<std::vector<double>> obstacles; 
    FG_eval(const Eigen::VectorXd& coeff, const vector<std::vector<double>>& obstacles){
        this->coeff = coeff;
        this->obstacles = obstacles; // in local frame(x0)
    };
    
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void setSpeedProfile(std::vector<double>& speed_profile){
        double x1 = 0.0, x2, y1, y2; 
        for(int i = 0; i < N; i++){
            y1 = coeff[0] + coeff[1]*x1 + coeff[2]*pow(x1,2) + coeff[3]*pow(x1,3);
            x2 = x1 + _dx; 
            y2 = coeff[0] + coeff[1]*x2 + coeff[2]*pow(x2,2) + coeff[3]*pow(x2,3);
            double dx = x2-x1, dy = y2-y1; 
            double ang = atan2(dy, dx); 
            if(dx != 0 && dy != 0){
                double dangle = abs(pi2pi(ang-(coeff[1] + 2*coeff[2]*x1 + 3*coeff[3]*pow(x1,2))));
                if(dangle >= 0.4363) speed_profile[i] = -target_speed;
            }
        }
        speed_profile.back() = 0.0; 
    }

    void operator()(ADvector& fg, const ADvector& var){

        // fg[0] = 0;
        std::vector<double> speed_profile(N,target_speed);
        setSpeedProfile(speed_profile); 
        for(int t = 0;t<N;t++){

            //cost for cte
            fg[0] += 1000*CppAD::pow(var[cte_start+t],2);
            // cost for orientation errror
            fg[0] += 20000*CppAD::pow(var[we_start+t],2);
            //cost of target speed
            // fg[0] += 10*CppAD::pow((var[v_start+t]-speed_profile[t]),2);
            fg[0] += 10*CppAD::pow((var[v_start+t]-target_speed),2);
            // cost of obstacles
            // for(int ii = 0; ii < obstacles.size(); ii++)
            //     fg[0] += 10000*CppAD::pow(1/sqrt(pow(var[x_start+t]-obstacles[ii][0],2)+pow(var[y_start+t]-obstacles[ii][1],2)),2);
            // TODO add actuator cost

            //TODO  add shock cost
            // dilute actuators
            // if (t < N - 1) {
            //     // steering
            //     fg[0] += CppAD::pow(var[delta_start + t], 2);
            //     // throttle/brake
            //     fg[0] += CppAD::pow(var[a_start + t], 2);
            // }
            
            // // actuator shock absorber
            // if (t < N - 2) {
            //     //steering
            //     fg[0] += 30000.0*CppAD::pow(var[delta_start + t + 1] - var[delta_start + t], 2);
            //     //fg[0] += CppAD::pow(var[delta_start + t + 1] - var[delta_start + t], 2);

            //     //acceleration
            //     fg[0] += CppAD::pow(var[a_start + t + 1] - var[a_start + t], 2);    
            // }

        }

        // constraint for the initial state
        fg[1+x_start] = var[x_start];
        fg[1+y_start] = var[y_start];
        fg[1+w_start] = var[w_start];
        fg[1+v_start] = var[v_start];
        fg[1+cte_start] = var[cte_start];
        fg[1+we_start] = var[we_start];
        
        
        for(int ii = 0; ii < obstacles.size(); ii++)
            fg[1+obs_start+ii*N] = sqrt(pow(var[x_start]-obstacles[ii][0],2)+pow(var[y_start]-obstacles[ii][1],2));
        

        // now we make sure that x-->x+1 can be reached through the integration of x 
        // this part has the dynamics and integration combined
        for(int t=1;t<N;t++){
            //current state
            AD<double> x0 = var[x_start+t-1];
            AD<double> y0 = var[y_start+t-1];
            AD<double> w0 = var[w_start+t-1];
            AD<double> v0 = var[v_start+t-1];
            AD<double> cte0 = var[cte_start+t-1];
            AD<double> we0 = var[we_start+t-1];
            AD<double> delta0 = var[delta_start+t-1];
            AD<double> a0 = var[a_start+t-1];
            //next state
            AD<double> x1 = var[x_start+t];
            AD<double> y1 = var[y_start+t];
            AD<double> w1 = var[w_start+t];
            AD<double> v1 = var[v_start+t];
            AD<double> cte1 = var[cte_start+t];
            AD<double> we1 = var[we_start+t];
            // AD<double> delta1 = var[delta_start+t];
            // AD<double> a1 = var[a_start+t];

            //calc ref point
            AD<double> f0 = coeff[0] + coeff[1]*x0 + coeff[2]*pow(x0,2) + coeff[3]*pow(x0,3);
            //cacl ref point angle
            AD<double> w_ref = CppAD::atan(coeff[1] + 2*coeff[2]*x0 + 3*coeff[3]*pow(x0,2));
            // now we say that the next state - integrated current state should be 0
            fg[1+x_start+t] = x1 - (x0 +v0*CppAD::cos(w0)*dt);
            fg[1+y_start+t] = y1 - (y0 +v0*CppAD::sin(w0)*dt);
            fg[1+w_start+t] = w1 - (w0 +v0*delta0/Lf*dt);
            fg[1+v_start+t] = v1 - (v0 + a0*dt);
            fg[1+cte_start+t] = cte1 - ((f0-y0)+(v0*CppAD::sin(we0)*dt));
            fg[1+we_start+t] = we1 - ((w0-w_ref)+v0*delta0/Lf*dt);
            for(int ii = 0; ii < obstacles.size(); ii++){
                // cout << sqrt(CppAD::pow(x1-obstacles[ii][0],2)+CppAD::pow(y1-obstacles[ii][1],2)) << endl;
                AD<double> obsx = obstacles[ii][0];
                AD<double> obsy = obstacles[ii][1];
                fg[1+obs_start+(ii*N)+t] = CppAD::sqrt(CppAD::pow(x1-obsx,2)+CppAD::pow(y1-obsy,2));
                // if(fg[1+obs_start+(ii*N)+t] < 2.0) cout << t << " ERROROROROROROROR" << endl;
            }
        }

        //std::cout<<"inside FG_val";

    }
};


vector<double> MPC::mpc_solve(Eigen::VectorXd state,Eigen::VectorXd coeff,
                              std::vector<std::vector<double>>& obstacles){
    size_t i;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    size_t n_var = (N)*(6)+(N-1)*2; // number of states for all time + number of control i/p for time-1
    size_t n_constraints = (N)*(6+obstacles.size()); //to ensure that integration contraints are meet
    Dvector var(n_var);
    for(int i = 0;i<n_var;i++){
        var[i] = 0;
    }
    
    //copying states into CPPAD vector
    var[x_start] = state[0];
    var[y_start] = state[1];
    var[w_start] = state[2];
    var[v_start] = state[3];
    var[cte_start] = state[4];
    var[we_start] = state[5];
    //upper and lower limits for var(states)
    Dvector var_lowerbounds(n_var);
    Dvector var_upperbounds(n_var);

    for (int i = 0;i<delta_start;i++){
        var_lowerbounds[i] = -99999999.;
        var_upperbounds[i] = 99999999.;
    }

    for (int i = v_start; i < cte_start; i++){
        var_lowerbounds[i] = 0.5;
    }

    for(int i = delta_start;i<a_start;i++){
        var_lowerbounds[i] = -0.4363*1.5;
        var_upperbounds[i] =  0.4363*1.5; // 25 deg to rad
    }

    for(int i=a_start;i<n_var;i++){
        var_lowerbounds[i] = -1.;
        var_upperbounds[i] = 1.;
    }

    Dvector constraints_lb(n_constraints); 
    Dvector constraints_up(n_constraints);
    for(int i=0;i<n_constraints;i++){
        constraints_up[i] = 0.;
        constraints_lb[i] = 0.;
    }
    for(int i = N*(6); i < n_constraints; i++){
        constraints_up[i] = 1000.; 
        constraints_lb[i] = 1.;   
    }

    constraints_lb[x_start] = state[0];
    constraints_lb[y_start] = state[1];
    constraints_lb[w_start] = state[2];
    constraints_lb[v_start] = state[3];
    constraints_lb[cte_start] = state[4];
    constraints_lb[we_start] = state[5];

    constraints_up[x_start] = state[0];
    constraints_up[y_start] = state[1];
    constraints_up[w_start] = state[2];
    constraints_up[v_start] = state[3];
    constraints_up[cte_start] = state[4];
    constraints_up[we_start] = state[5];
    
    FG_eval fg_eval(coeff, obstacles);
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          2\n";

    // std::string options;

    options += "Sparse true             forward\n";
    options += "Sparse true             reverse\n";
    options += "Numeric max_cpu_time    0.5\n";

    CppAD::ipopt::solve_result<Dvector> result;
    // std::cout<<"calling ipopt";
    //********* calling the ipopt solver **********//
    CppAD::ipopt::solve(options,var,var_lowerbounds,var_upperbounds,constraints_lb,constraints_up,fg_eval,result);

    // // Parsing the result
    std::cout << "Solution Status :" << result.status << "\n";
    // std::cout << "cost :" << result.obj_value << "/n";
    // return control at the 1st timestamp

    bool ok = true;
    ok &= result.status == CppAD::ipopt::solve_result<Dvector>::success;
    std::vector<double> control;
    control.push_back(result.x[delta_start]);
    control.push_back(result.x[a_start]);
    for(int i = 0; i < N; i++){
        control.push_back(result.x[x_start+i]); 
        control.push_back(result.x[y_start+i]);
    }
    std::cout<<"delta " << control[0] << std::endl;
    std::cout<<"a " << control[1] << std::endl;
    auto cost = result.obj_value; 
    // cout << "cost: " << cost << endl;
    return control;

}





// int main (){
//     // generate trajectory
//     std::vector<double> cx((60.0)/_dx,0.0);
//     std::vector<double> cy(cx.size(), 0.);
//     std::vector<double> cd(cx.size(), 0.);
//     for(int i = 1; i < cx.size();i++){
//         cx[i] = cx[i-1]+0.1;
//         cy[i] = sin(cx[i]/5.0);//*cx[i]/2.0;//sin(ix / 5.0) * ix / 2.0 for ix in cx]
//         cd[i] = sqrt(0.01 + pow(cy[i]-cy[i-1], 2)); 
//     }
    
//     // obstacles
//     std::vector<std::vector<double>> obstacles, obstacles_local;
//     // std::vector<double> tmp  = {16.0, 0.0}, tmp2 = {};
//     // obstacles.push_back(tmp); 
//     obstacles = {{16.0, 0.0}, {40.0, 1.2}, {26.0, 0.0}};
//     obstacles_local = obstacles;

//     //initalize state
//     // double _x, _y; 
//     // std::cout << "enter x: " << endl;
//     // std::cin >> _x;
//     // std::cout << "enter y: " << endl;
//     // std::cin >> _y; 
//     state veh_ (0,3,0,2.0);
//     // set speed
//     std::vector<state> veh_state;
//     veh_state.push_back(veh_);
//     std::vector<double> x,y,w,v,t;

//     int max_time = 100;
//     double time = 0, delta = 0, a =1;
//     double lookaheadDist = 0; 
//     //std::cout << "before while ";

//     MPC controller;
//     int fit_numberof_pts = 14;
//     Eigen::VectorXd x_14pts(fit_numberof_pts);
//     Eigen::VectorXd y_14pts(fit_numberof_pts);
//     int last_target_idx = 0;
//     while (time < max_time ){
//         transform_to_local(veh_, obstacles, obstacles_local);
//         int clst_indx = calc_lookahead_pt(cx, cy,veh_); // find the closest point

//         get_14points(clst_indx, x_14pts, y_14pts, cx, cy, cd, lookaheadDist); //find the closest 14 points
//         std::vector<double> x14pts(x_14pts.size(),0),y14pts(x_14pts.size(),0);
//         for(int i=0;i<x_14pts.size();i++){
//             x14pts[i] = x_14pts[i];
//             y14pts[i] = y_14pts[i];;

//         }
//         transform_to_local(veh_,x_14pts,y_14pts);
//         //std::cout<<"14 points";
//         Eigen::VectorXd coeff;
//         try{
//             polyfit(coeff, x_14pts,y_14pts,3);
//         }catch(int e){
//             cout << "DONE" << endl;
//             break;
//         }
//         vector<double> xplot(1,x_14pts[0]);
//         vector<double> yplot;
//         yplot = polyvec(coeff, xplot, 80);
//         transform_to_global(veh_, xplot, yplot);

//         // std::cout<<"start coeff";
//         // for(int i=0;i<coeff.size();i++){
//         //     std::cout<<coeff[i]<<std::endl;
//         // }
//         // std::cout<<"start coeff";

//         double cte = coeff[0];
//         // std::cout<<"CTE"<<cte;
//         int rand;
//         // std::cin>> rand;
//         double we = -atan(coeff[1]);
//         // std::cout<<we;
//         Eigen::VectorXd vehstate(6);
//         vehstate << 0, 0, 0, veh_.v, cte, we;
        
//         auto res = controller.mpc_solve(vehstate,coeff, obstacles_local);
//         veh_.update(res[1],res[0],dt);
//         cout << "*********" << endl;
//         for(auto i : obstacles){            
//             cout << sqrt(pow(veh_.x-i[0],2) + pow(veh_.y-i[1],2)) << endl;
//         }
//         cout << "*********" << endl << endl;
//         vector<double> xmpc, ympc; 
//         // cout << "size " << res.size() << endl;
        
//         double xtmp = 0, ytmp = 0; 
//         lookaheadDist = 0; 
//         for(int i = 2; i < res.size(); i += 2){
//             lookaheadDist += sqrt(pow(res[i]-xtmp,2) + pow(res[i+1]-ytmp,2)); 
//             xtmp = res[i]; 
//             ytmp = res[i+1]; 
//             xmpc.push_back(xtmp); 
//             ympc.push_back(ytmp); 
//         }
//         transform_to_global(veh_, xmpc, ympc);

//         time = time + dt;
//         veh_state.push_back(veh_);
//         last_target_idx++;
//         if(last_target_idx == cx.size()-14)
//             break;
//         #if PLOT
//             x.push_back(veh_state.back().x);
//             y.push_back(veh_state.back().y);
//             // Clear previous plot
// 			plt::clf();
// 			// Plot line from given x and y data. Color is selected automatically.            
// 			plt::named_plot("Car",x,y,"*k");
//             plt::named_plot("Track",cx,cy,"-r");   
//             plt::named_plot("poly_points",x14pts,y14pts,"-b");
//             /****/
//             plt::named_plot("mpc_traj", xmpc, ympc, "go");
//             // plt::named_plot("poly", xplot, yplot, "r+");
//             for(int i = 0; i < obstacles.size(); i++) plotCircle(obstacles[i][0], obstacles[i][1], 1.0);
//             /****/
//             plt::legend();
//             plt::axis("equal");
//             // plt::xlim(-1,60);
//             // plt::plot(cx,cdd,"-b");
//             plt::pause(0.1);

//         #endif
//     }
//     plt::show(); 
//     return 0;
// }

//g++ mpc_steering.cpp -I/usr/include/python2.7 -lpython2.7 -lipopt -std=c++11
