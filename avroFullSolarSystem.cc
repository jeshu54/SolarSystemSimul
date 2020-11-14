// Authors: Brian Ng and Avro Mukherjee

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <random>
#include <cmath>
using namespace std;

const double G = 6.67408E-11; // Gravitational constant
constexpr double pi = 3.14159265358979;
const double M = 1.9891e30; // Mass of the sun

default_random_engine generator(0); // seed with 0 for repeatable pseudorandom values
//default_random_engine generator; // replace the previous with this line for different results every time

class Vec3d {
private:
  double x, y, z;
public:
  Vec3d(double x, double y, double z) : x(x), y(y), z(z) {}

  double getX() {
    return x;
  }
  
  double getY() {
    return y;
  }

  double getZ() {
    return z;
  }

  void setX(double x) {
    this->x = x;
  }

  void setY(double y) {
    this->y = y;
  }

  void setZ(double z) {
    this->z = z;
  }

  friend ostream& operator <<(ostream& s, const Vec3d& v) {
    return s << '(' << v.x << ", " << v.y << ", "  << v.z << ")";
  }
};

class Body {
private:
  string name;
  double mass;
  Vec3d pos; // positions are random
  Vec3d v;
  Vec3d a;
public:
  Body(string name, double mass, Vec3d pos, Vec3d v, Vec3d a) : name(name), mass(mass), pos(pos), v(v), a(a) {}

  string getName() {
    return name;
  }

  double getMass() {
    return mass;
  }
  
  Vec3d getPos() {
    return pos;
  }

  Vec3d getV() {
    return v;
  }

  Vec3d getA() {
    return a;
  }

  friend ostream& operator <<(ostream& s, const Body& b) {
    return s << b.name << ' ' << b.pos << ' ' << b.v << ' ' << b.a;
  }
};

class SolarSystem {
private:
  string fileName;
  string line, name, orbits;
  double mass, diam, perihelion, aphelion, avgDist;
  double xCoord, yCoord, zCoord, vel, velX, velY, velZ, acc, accX, accY, accZ;
  double theta = pi;
  double positionX [4] = {0,-5.88841e+10,-1.08593e+11,-1.48706e+11}; // initial value will be taken from bodies[i].getPos().getX() for each object, we can hardcode it. Also for storing new values.
  
  double positionY [4] = {0,0.00019026,0.000350875,0.000480481}; // initial value will be taken from bodies[i].getPos().getY() for each object, we can hardcode it. Also for storing new values.

  vector<Body> bodies;  // all bodies in the solar system
public:
  SolarSystem(string fileName) : fileName(fileName) {
    ifstream f(fileName);

    if(!f.is_open()) {
  	  cout << "The file is not open." << endl;
    }
    else {
      while (getline(f, line)) {
        if (line[0] == 'N' && line[1] == 'a' && line[2] == 'm') 
    	    continue;
        else {
          if (line[0] == 'M' && line[1] == 'o' && line[2] == 'o') {
            break;
          }
          else {
          stringstream sso(line);
          
          sso >> name >> orbits >> mass >> diam >> perihelion >> aphelion;
          avgDist = ((perihelion + aphelion)/2);

          randPosition();
          initVel();
          initAcc();

          Vec3d pos(xCoord, yCoord, zCoord);
          Vec3d v(velX, velY, velZ);
          Vec3d a(accX, accY, accZ);

          Body b(name, mass, pos, v, a);

          bodies.push_back(b);
          }
        }
      }
    }

    cout << "Initial Location\n";
    for(int i = 0; i < 4; i++) {
      cout << bodies[i].getName() << " (" << bodies[i].getPos().getX() << ", " << bodies[i].getPos().getY() << ", " << bodies[i].getPos().getZ() << ")\n";
    }
  }

  double randPosition() {
    uniform_real_distribution<double> randomDist(aphelion, perihelion);
    double r = randomDist(generator);

    if (name != "Sun") {
      xCoord = r * cos(theta);
      yCoord = r * sin(theta);
      zCoord = 0;
    }
    else {
      xCoord = 0;
      yCoord = 0;
      zCoord = 0;
    }

    return xCoord, yCoord, zCoord;
  }

  double initVel() {
    if (perihelion != 0) {
      vel = (sqrt((G*M)/(diam/2)))/100;
      //cout << vel << '\n';

      velX = vel * cos(pi);
      velY = vel * sin(pi);
      velZ = 0;
    }
    else {    // Sun is 0 everything
      velX = 0;
      velY = 0;
      velZ = 0;
    }
    return velX, velY, velZ;
  }

  double initAcc() {
    // acc = pow(vel, 2)/(diam/2);
    if (perihelion != 0) {
      acc = (pow(sqrt((G*M)/(diam/2)), 2)/(diam/2))/100;
      //cout << acc << '\n';

      accX = acc * cos(pi);
      accY = acc * sin(pi);
      accZ = 0;
    }
    else {    // Sun is 0 everything
      accX = 0;
      accY = 0;
      accZ = 0;
    }
    
    return accX, accY, accZ;
  }

  double accGrav(const double G, double m2, double x1, double y1, double z1, double x2, double y2, double z2) { // tells us how much m2 will accelerate us, m2 = different bodies in the system
    double dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
    double rsq = dx*dx + dy*dy + dz*dz;
    return G*m2/rsq;
  }
  
  void timeStep(double dt) {
    accX = 0;
    accY = 0;
    accZ = 0;
    theta += 2*pi/dt;

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if(bodies[i].getName() != bodies[j].getName()) {
          double acc = accGrav(G, bodies[j].getMass(), bodies[i].getPos().getX(), bodies[i].getPos().getY(), bodies[i].getPos().getZ(), bodies[j].getPos().getX(), bodies[j].getPos().getY(), bodies[j].getPos().getZ());
        //double acc = accGrav(G, bodies[j].getMass(), positionX[i], positionY[i], bodies[i].getPos().getZ(), positionX[j], positionY[j], bodies[j].getPos().getZ());
          accX += acc * cos(theta);
          accY += acc * sin(theta);
        }
        else {
          continue;
        }
      }

      if (bodies[i].getName() != "Sun") {
        //positionX[i] += (positionX[i] + bodies[i].getV().getX() * (0.5) * dt + accX * dt*dt); // To print position X for body i after acceleration due to all j bodies are calcualated
        //positionY[i] += (positionY[i] + bodies[i].getV().getY() * (0.5) * dt + accY * dt*dt); // To print position Y for body i after acceleration due to all j bodies are calcualated
        positionX[i] += (bodies[i].getPos().getX() + bodies[i].getV().getX() * (0.5) * dt + accX * dt*dt);
        positionY[i] += (bodies[i].getPos().getY() + bodies[i].getV().getY() * (0.5) * dt + accY * dt*dt);
      }
      else {
        positionX[i] = 0;
        positionY[i] = 0;
        continue;
      }
    accX = 0;
    accY = 0;
    }
  }

#if 0
  friend ostream& operator <<(ostream& s, const SolarSystem& ssystem) {
    for(int i = 0; i < 4; i++) {
      s << bodies[i].getName() << " (" << positionX[i] << ", " << positionY[i] << ", 0)\n";
    }
    return s;
  }
#endif

  void print() {
    cout << "Final Positions\n";
    for(int i = 0; i < 4; i++) {
      cout << bodies[i].getName() << " (" << positionX[i] << ", " << positionY[i] << ", 0)\n";
    }
  }
};

int main() {
  SolarSystem s("solarsystem.dat"); // this must not change
  // each body should have a random position in a circular or elliptical orbit around the sun
  // each body should start with an appropriate velocity
  // each body should have velocity = https://en.wikipedia.org/wiki/Circular_orbit v = sqrt(GM/r)

  double earthYear = 365.2425 * 24 * 60 * 60;
  const int numTimeSteps = 1000;
  double dt = earthYear / numTimeSteps; // dt = 31556.952

  for (int i = 0; i < numTimeSteps; i++) {
    s.timeStep(dt);
  }
  s.print();
  //cout << s; // print out the state of the solar system (where each object is, location only)
  // if you do it right, earth should be in roughly the same place...
}