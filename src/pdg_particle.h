#ifndef PDG_PARTICLE_H
#define PDG_PARTICLE_H

#pragma once
// A modified version of Andrea's code
#include <string>
#include <filesystem>
#include "interfaces.h"
#include "utils.h"
namespace powerhouse
{
    const std::string DATABASE = "pdg_database/baryons_mesons.dat";

    class pdg_particle : public I_particle
    {
    private:
        int _id;
        double _mass;
        std::string _name;
        float _spin;
        double _q;
        double _b;
        double _s;
        bool _particle;

    public:
        pdg_particle() : _id(0) {}
        pdg_particle(int id)
        {
            if (!std::filesystem::exists(DATABASE))
            {
                throw std::runtime_error("PDG database was not found.");
            }

            std::ifstream input(DATABASE);
            if (!input.is_open())
            {
                throw std::runtime_error("Failed to open PDG database.");
            }
            pdg_particle temp;
            _particle = id > 0;
            id = abs(id);
            std::string line;
            while (std::getline(input, line))
            {
                if (line.empty() || line[0] == '#')
                {
                    continue; // skip empty or comment lines
                }
                std::istringstream iss(line);
                if (!iss.fail())
                {
                    iss >> temp;

                    if (temp._id == id)
                    {
                        this->_b = temp._b * (_particle ? 1 : -1);
                        ;
                        this->_id = id;
                        this->_mass = temp._mass;
                        this->_name = temp._name;
                        this->_q = temp._q * (_particle ? 1 : -1);
                        this->_s = temp._s * (_particle ? 1 : -1);
                        this->_spin = temp._spin;
                        if (!_particle)
                        {
                            _name = _name.insert(0, "anti-");
                        }
                        break;
                    }
                }
            }
            if (temp._id == 0)
            {
                throw std::runtime_error("Particle not found");
            }
        }

        pdg_particle(pdg_particle &other)
        {
            this->_b = other._b;
            this->_id = other._id;
            this->_mass = other._mass;
            this->_name = other._name;
            this->_q = other._q;
            this->_s = other._s;
            this->_spin = other._spin;
        };
        ~pdg_particle() override
        {
        }
        std::string name() override { return _name; };
        double mass() override { return _mass; };
        int pdg_id() override { return _id; };
        double Q() override { return _q; };
        double B() override { return _b; };
        double S() override { return _s; };
        float spin() override { return _spin; };
        bool isparticle() override { return _particle; };

        friend std::istream &operator>>(std::istream &stream, pdg_particle &particle)
        {
            stream >> particle._id >> particle._mass >> particle._name >> particle._q >> particle._spin >> particle._b >> particle._s;
            return stream;
        };
        friend std::ostream &operator<<(std::ostream &stream, const pdg_particle &particle)
        {
            stream << particle._name << ", M = " << particle._mass
                   << ", (B,Q,S) =  (" << particle._b << "," << particle._q << "," << particle._s
                   << ") Spin-" << particle._spin;
            return stream;
        };
    };
}

#endif