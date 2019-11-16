/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <valarray>

#include "helper_functions.h"
#include "multiv_gauss.h" 

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 200;  // TODO: Set the number of particles

  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(x,std[0]);
  std::normal_distribution<double> dist_y(y,std[1]);
  std::normal_distribution<double> dist_theta(theta,std[2]);

  for (int i = 0; i< num_particles; ++i) {
    double sample_x, sample_y, sample_theta;

    sample_x = dist_x(generator); 
    sample_y = dist_y(generator);
    sample_theta = dist_theta(generator);

    Particle new_part;
    new_part.id = i;
    new_part.x = sample_x;
    new_part.y = sample_y;
    new_part.theta = sample_theta;
    new_part.weight = 1.0;
    
    particles.push_back(new_part);
    weights.push_back(1.0);

  is_initialized = true;

  //std::cout << "Initialised" << std::endl;

  return;
  } 

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine pred_gen;
  //std::cout << "predict" << std::endl;
  
  for (std::size_t i=0; i<particles.size(); ++i) {
    double update_x, update_y, update_theta;

    double p_x, p_y, p_theta;
    p_x = particles[i].x;
    p_y = particles[i].y;
    p_theta = particles[i].theta;
    
    if (abs(yaw_rate) > 1e-5) {
      update_theta = p_theta + yaw_rate*delta_t;
      update_x = p_x + velocity/yaw_rate * (std::sin(p_theta + yaw_rate*delta_t) - std::sin(p_theta) ); 
      update_y = p_y + velocity/yaw_rate * (std::cos(p_theta) - std::cos(p_theta + yaw_rate*delta_t) ); 
      
    } else {
      update_x = p_x + velocity * delta_t * std::cos(p_theta); 
      update_y = p_y + velocity * delta_t * std::sin(p_theta);
      update_theta = p_theta;
    }
    
    std::normal_distribution<double> dist_x(update_x,std_pos[0]);
    std::normal_distribution<double> dist_y(update_y,std_pos[1]);
    std::normal_distribution<double> dist_theta(update_theta,std_pos[2]);

    particles[i].x = dist_x(pred_gen);
    particles[i].y = dist_y(pred_gen);
    particles[i].theta = dist_theta(pred_gen);

  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // predicted the landmark list | observations is the transformed observations


//vector<LandmarkObs> associate_obs(vector<LandmarkObs> predicted, 
//                                     vector<LandmarkObs>& observations) {
                                  
  std::cout << "running associations" << std::endl;

  for (std::size_t i=0; i<observations.size(); ++i) {
    
    LandmarkObs cur_obs = observations.at(i);

    double min_dist=std::numeric_limits<double>::max();

    for (std::size_t j = 0; j<predicted.size(); ++j) {
      LandmarkObs cur_pred = predicted[j];

      double error_dist = dist(cur_pred.x, cur_pred.y, cur_obs.x, cur_obs.y);
      
      if (error_dist < min_dist) {
        min_dist = error_dist;
        observations.at(i).id = j;
      }
    }

  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double sig_x, sig_y;
  sig_x = std_landmark[0];
  sig_y = std_landmark[1];

  for (std::size_t i = 0; i < particles.size(); ++i) {
    
    Particle cur_part = particles.at(i);

    vector<LandmarkObs> trans_obs;
    double weight_update = 1.0;

    int obs_length = observations.size();

    for (int j = 0; j < obs_length; ++j) {
      LandmarkObs cur_obs = observations.at(j);

      LandmarkObs trans_observation;
      trans_observation.id = cur_obs.id;
      trans_observation.x = cur_part.x + (cos(cur_part.theta) * cur_obs.x) - (sin(cur_part.theta) * cur_obs.y);
      trans_observation.y = cur_part.y + (sin(cur_part.theta) * cur_obs.x) + (cos(cur_part.theta) * cur_obs.y);

      trans_obs.push_back(trans_observation);
    }
    
    std::vector<LandmarkObs> close_landmarks_lst;

    std::cout << "observations list: " << observations.size() << std::endl;
    std::cout << "transformed observations list: " << trans_obs.size() << std::endl;
    
    // find observations within radius of sensor

    for (std::size_t map_it = 0; map_it < map_landmarks.landmark_list.size(); ++map_it) {
      
      Map::single_landmark_s cur_landmark = map_landmarks.landmark_list[map_it];

      double distance = dist(cur_part.x, cur_part.y, cur_landmark.x_f, cur_landmark.y_f);
      if (distance < sensor_range) {
        LandmarkObs cur_lmk;
        cur_lmk.id = cur_landmark.id_i;
        cur_lmk.x = cur_landmark.x_f;
        cur_lmk.y = cur_landmark.y_f;

        // just being lazy don't wanna edit dataAssociation again
        close_landmarks_lst.push_back(cur_lmk);
      } 

    }

    std::cout << "close landmark's list: " << close_landmarks_lst.size() << std::endl;
    
    // this bit breaks
    // we aren't calculating weight properly
    dataAssociation(close_landmarks_lst, trans_obs);
    // after this step all the trans_obs will have a landmark id and a 

    double weight_prob;
    for (std::size_t k; k < trans_obs.size(); ++ k) {

      LandmarkObs cur_obs = trans_obs[k];
      LandmarkObs close_landmark = close_landmarks_lst[cur_obs.id];

      weight_prob = multiv_prob(sig_x, sig_y, cur_obs.x, cur_obs.y, close_landmark.x, close_landmark.y);
      std::cout << "weight_prob "  << weight_prob << std::endl;

      weight_update *= weight_prob;
      
    }

    std::cout << "weight_update " << weight_update << std::endl;
    // calculate the new weight
    particles[i].weight = weight_update;
    weights[i] = particles[i].weight;
  
  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::cout << "Resample" << std::endl;
  
  std::default_random_engine generator;

  std::vector<Particle> result;

  for (std::size_t i = 0; i < weights.size(); ++i ) {

    std::discrete_distribution<int> distrib(weights.begin(), weights.end());
    int value = distrib(generator);
    result.push_back(particles[value]);
    //result.at(i) = particles.at(value);
  }

  particles = result;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}