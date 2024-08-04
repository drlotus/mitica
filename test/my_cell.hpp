class my_cell : public hydro::Icell<utils::geometry::four_vector, utils::r2_tensor>
{
public:
    my_cell(){};
    my_cell(const my_cell &other)
    {
        _tau = other._tau;
        _x = other._x;
        _y = other._y;
        _eta = other._eta;
        _dsigma = other._dsigma;
        _u = other._u;
        _T = other._T;
        _mub = other._mub;
        _muq = other._muq;
        _mus = other._mus;
        _dbeta = other._dbeta;
        _du = other._du;
    }
    ug::four_vector milne_coords() const override { return ug::four_vector({_tau, _x, _y, _eta}, false); }
    ug::four_vector thermodynamics() const override { return ug::four_vector({_T, _mub, _muq, _mus}, false); }
    utils::r2_tensor du_ll() const override { return _du; }
    utils::r2_tensor dbeta_ll() const override { return _dbeta; }
    ug::four_vector four_vel() const override { return _u; }
    ug::four_vector disgma() const override { return _dsigma; }
    bool is_spacelike() override { return normal_sq() < 0; }
    ug::four_vector acceleration() override
    {
        if (!_acc)
        {
            calculte_ac();
        }
        return *_acc;
    }
    double u_dot_n() override
    {
        return _u * _dsigma;
    }

    friend std::istream &operator>>(std::istream &stream, my_cell &cell)
    {
        stream >> cell._tau >> cell._x >> cell._y >> cell._eta;
        stream >> cell._dsigma;
        stream >> cell._u;
        stream >> cell._T >> cell._mub >> cell._muq >> cell._mus;

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                stream >> cell._dbeta[i][j];
            }
        }
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                stream >> cell._du[i][j];
            }
        }

        return stream;
    }
    friend std::ostream &operator<<(std::ostream &stream, const my_cell &cell)
    {
        stream << cell.milne_coords() << "\t";
        stream << "T = " << cell._T << "\t mu_B = " << cell._mub << "\t mu_Q = " << cell._muq << "\t mu_S = " << cell._mus;
        return stream;
    }

    std::ostream &write(std::ostream &os) override
    {
        os << *this;
        return os;
    }

    double normal_sq() override
    {
        if (!_normal_size)
        {
            _normal_size = std::make_unique<double>(_dsigma.norm_sq());
        }
        return *_normal_size;
    }

private:
    double _tau, _x, _y, _eta;
    ug::four_vector _u;
    ug::four_vector _dsigma;
    double _T, _mub, _muq, _mus;
    utils::r2_tensor _dbeta;
    utils::r2_tensor _du; // derivatives of the 4-velocity in Cartesian coordinates
    std::unique_ptr<double> _normal_size;
    std::unique_ptr<double> _theta;
    std::unique_ptr<double> _b_theta;
    std::unique_ptr<utils::r2_tensor> _delta_ll;
    std::unique_ptr<utils::r2_tensor> _delta_ul;
    std::unique_ptr<utils::r2_tensor> _delta_uu;
    std::unique_ptr<ug::four_vector> _acc;
    std::unique_ptr<utils::r2_tensor> _th_vorticity;
    std::unique_ptr<utils::r2_tensor> _th_shear;
    std::unique_ptr<utils::r2_tensor> _shear;
    std::unique_ptr<ug::four_vector> _f_vorticity_vec;
    std::unique_ptr<utils::r2_tensor> _f_vorticity;
    std::unique_ptr<utils::r2_tensor> _gradu;
    void calculate_shear() {}
    void calculate_fvorticity_vec() {}
    void calculate_fvorticity() {}
    void calculate_th_vorticity() {}
    void calculate_th_shear() {}
    void calculte_ac()
    {

        ug::four_vector _(false);
        for (size_t i = 0; i < 4; i++)
        {
            _[i] = 0;
            for (size_t j = 0; j < 4; j++)
            {
                _[i] += _u[j] * _du[j][i] * utils::gmumu[i];
            }
        }
        _acc = std::make_unique<ug::four_vector>(_);
    }

    std::unique_ptr<double> _sigma_norm;
    std::unique_ptr<double> _fvort_norm;
    std::unique_ptr<double> _tvort_norm;
    std::unique_ptr<double> _tshear_norm;
    std::unique_ptr<double> _acc_norm;
};