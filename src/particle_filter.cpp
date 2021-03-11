/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 * Updated on: Mar 9, 2021 by Albert Fu
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

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
 // std::cout << "Initialization Started" << std::endl;
  num_particles = 5;  // TODO: Set the number of particles
  // Creates a normal (Gaussian) distribution for x, y, theta
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle temp;
    temp.id = i;
    temp.x = dist_x(gen);
    temp.y = dist_y(gen);
    temp.theta = dist_theta(gen);
    temp.weight = 1;
    particles.push_back(temp);
  }
  is_initialized = true;
 // std::cout << "Initialization Succeed" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    //std::cout << "Prediction Started" << std::endl;
    std::default_random_engine gen;
	for (int i = 0; i < particles.size(); i++) {
      if(yaw_rate > 0.000001){
      normal_distribution<double> dist_x(particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+delta_t*yaw_rate)-sin(particles[i].theta)), std_pos[0]);
      normal_distribution<double> dist_y(particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+delta_t*yaw_rate)), std_pos[1]);
      normal_distribution<double> dist_theta(particles[i].theta+delta_t*yaw_rate, std_pos[2]);
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
      }
      else
      {
      normal_distribution<double> dist_x(particles[i].x + velocity*cos(particles[i].theta)*delta_t, std_pos[0]);
      normal_distribution<double> dist_y(particles[i].y + velocity*sin(particles[i].theta)*delta_t, std_pos[1]);
      normal_distribution<double> dist_theta(particles[i].theta+delta_t*yaw_rate, std_pos[2]);
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
      }
  	}
 //   std::cout << "Prediction Succeed" << std::endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations, double sensor_range) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
 // std::cout << "Data Association Started" << std::endl;
  vector<int> matched_mark;
  int res = 0;
  for (int i = 0; i < observations.size(); i++){
     res = sensor_range;
     observations[i].id = -1;
  	for(int j = 0; j < predicted.size(); j++){
      if(j >= observations.size())
        {
          break;
        }
        if (std::find(matched_mark.begin(), matched_mark.end(), j) != matched_mark.end())
      		{continue;}
      	if(dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y) < res){
            res = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
        	observations[i].id = j; 
        }  
      }
    matched_mark.push_back(observations[i].id);
  }
  /*std::cout << "Landmark:" << std::endl;
  for(int j = 0; j < predicted.size(); j++){
  std::cout << j << ": " << predicted[j].x << ", " << predicted[j].y << std::endl;
  }
  std::cout << "Observation:" << std::endl;
  for(int j = 0; j < observations.size(); j++){
  std::cout << j << ": " << observations[j].x << ", " << observations[j].y << ", " << observations[j].id << std::endl;
  }
  */
  //getchar();
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
   // std::cout << "Update Weights Started" << std::endl;
  	
  	vector<LandmarkObs> land_mark;
    LandmarkObs temp_land_mark;

 	double sum = 0;
    for(int i = 0; i < map_landmarks.landmark_list.size(); i++){
      temp_land_mark.x = map_landmarks.landmark_list[i].x_f;
      temp_land_mark.y = map_landmarks.landmark_list[i].y_f;
      temp_land_mark.id = map_landmarks.landmark_list[i].id_i;
      land_mark.push_back(temp_land_mark);
    }
  
	for (int i = 0; i < particles.size(); i++) {
      vector<LandmarkObs> temp;
      particles[i].weight = 1;
      vector<LandmarkObs> land_mark_in_range;
      for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
        if(dist(land_mark[j].x, land_mark[j].y, particles[i].x, particles[i].y) <= sensor_range){
        land_mark_in_range.push_back(land_mark[j]);
        }
      }
      for (int j = 0; j < observations.size(); j++){
        temp_land_mark.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
        temp_land_mark.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
        temp_land_mark.id = observations[j].id;
        temp.push_back(temp_land_mark);
      }
      dataAssociation(land_mark_in_range, temp, sensor_range);
      for (int j = 0; j < observations.size(); j++){
        if(j >= land_mark_in_range.size())
        {
          break;
        }
        /*if(multiv_prob(std_landmark[0], std_landmark[1], temp[j].x, temp[j].y, land_mark_in_range[temp[j].id].x, land_mark_in_range[temp[j].id].y) == 0)
        {
          std::cout << "temp[j].x: " << temp[j].x << std::endl;
          std::cout << "temp[j].y: " << temp[j].y << std::endl;
          std::cout << "land_mark_in_range[temp[j].id].x: " << land_mark_in_range[temp[j].id].x << std::endl; 
          std::cout << "land_mark_in_range[temp[j].id].y: " << land_mark_in_range[temp[j].id].y << std::endl; 
        }*/
        particles[i].weight *= multiv_prob(std_landmark[0], std_landmark[1], temp[j].x, temp[j].y, land_mark_in_range[temp[j].id].x, land_mark_in_range[temp[j].id].y);
      }
      sum += particles[i].weight;
  	}
  for (int i = 0; i < particles.size(); i++){
      if(sum == 0){
 /*      std::cout << "particles.weight: " << particles[0].weight << std::endl;
       std::cout << "land_mark_in_range.size(): " << land_mark_in_range.size() << std::endl;
       std::cout << "observations.size(): " << observations.size() << std::endl; 
       getchar(); */
      }
  	particles[i].weight = particles[i].weight/sum;
    //std::cout << "particles[i].weight: " << particles[i].weight << std::endl;
  }
 // std::cout << "Update Weights Succeed" << std::endl;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::uniform_real_distribution<> dis(0, 1.0);
//  std::cout << "particles.size(): " << particles.size() << std::endl;
  weights.clear();
 // std::cout << "Resample Started" << std::endl;
  for (int i = 0; i < particles.size(); i++) {
    weights.push_back(particles[i].weight);  
  }
 // std::cout << "particles.size(): " << particles.size() << std::endl;
  std::vector<Particle> temp;
  int index = round(dis(gen) * particles.size());
  //std::cout << "Current Index: " << index << std::endl;
  //getchar();
  double beta = 0;
  double mw = *max_element(weights.begin(), weights.end());
 // std::cout << "Start Resample Wheel" << std::endl;
  //std::cout << "particles.size(): " << particles.size() << std::endl;
  for (int i = 0; i < particles.size(); i++) {
   // std::cout << "Current Index: " << index << std::endl;
    beta += dis(gen) * 2 * mw;
    while (beta > weights[index]){
        beta -= weights[index];
        index = (index + 1) % particles.size();
    }
    temp.push_back(particles[index]);
  }
 // std::cout << "Resample Wheel Finished" << std::endl;
  particles = temp;
 /* for (int i = 0; i < particles.size(); i++) {
	particles[i].weight = 1.0;
  } */
   // std::cout << "Resample Succeed" << std::endl;
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

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}