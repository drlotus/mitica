#include <memory>
#include "vhlle_fcell.h"
#include "interfaces.h"
#include <map>
#pragma once
namespace hydro
{
    /// @brief singleton factory that is used for registeration and creation of analytical solution
    /// @tparam C the cell's type
    /// @tparam V the four-vector's type
    /// @tparam T the rank-2 tenosor's type
    template <typename C, typename V, typename T>
    class solution_factory
    {
    public:
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        using solution_creator = std::function<std::unique_ptr<I_solution<C, V, T>>()>;
        std::unique_ptr<I_solution<C, V, T>> create(const std::string &name)
        {
            const auto &it = _map.find(name);
            if (it != _map.end())
            {
                return it->second();
            }
            else
            {
                throw std::runtime_error("Solution not found!");
            }
        }
        void regsiter_solution(const std::string &name, solution_creator creator)
        {
            _map[name] = std::move(creator);
        }
        static std::shared_ptr<solution_factory<C, V, T>> &factory()
        {
            std::call_once(
                only_one,
                []()
                {
                    solution_factory<C, V, T>::_factory_instance.reset(new solution_factory());
                });
            return solution_factory<C, V, T>::_factory_instance;
        }
        ~solution_factory() { _map.clear(); }

    private:
        solution_factory() {}
        solution_factory(const solution_factory<C, V, T> &rs) = delete;
        solution_factory &operator=(const solution_factory<C, V, T> &rs) = delete;
        std::unordered_map<std::string, solution_creator> _map;
        static std::once_flag only_one;
        static std::shared_ptr<solution_factory<C, V, T>> _factory_instance;
    };
    template <typename C, typename V, typename T>
    std::shared_ptr<solution_factory<C, V, T>> solution_factory<C, V, T>::_factory_instance = nullptr;
    template <typename C, typename V, typename T>
    std::once_flag solution_factory<C, V, T>::only_one;
}

namespace powerhouse
{
    /// @brief singleton factory that is used for registeration and creation of calculators
    /// @tparam C the cell's type
    /// @tparam P particle's type
    /// @tparam O output's type
    template <typename C, typename P, typename O>
    class calculator_factory
    {
    public:
        static_assert(std::is_base_of<powerhouse::I_particle, P>::value, "P must inherit from I_particle");
        static_assert(is_template_base_of<hydro::I_cell, C>::value, "C must inherit from I_cell");
        static_assert(is_template_base_of<powerhouse::I_output, O>::value, "O must inherit from I_output");
        using calculator_creator = std::function<std::unique_ptr<I_calculator<C, P, O>>()>;
        std::unique_ptr<I_calculator<C, P, O>> create(const utils::program_options &options)
        {
            powerhouse::calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};
            
            const auto &it = _map.find(key);
            if (it != _map.end())
            {
                return it->second();
            }
            else
            {
                throw std::runtime_error("Calculator not found!");
            }
        }
        ~calculator_factory() { _map.clear(); }

        void register_calculator(utils::program_options options, calculator_creator creator)
        {
            powerhouse ::calculator_key key{options.program_mode, options.polarization_mode, options.yield_mode};
            _map[key] = std::move(creator);
        }

        static std::shared_ptr<calculator_factory<C, P, O>> &factory()
        {
            std::call_once(
                only_one,
                []()
                {
                    calculator_factory<C, P, O>::_factory_instance.reset(new calculator_factory<C, P, O>());
                });
            return calculator_factory<C, P, O>::_factory_instance;
        }

    private:
        calculator_factory() {}
        calculator_factory(const calculator_factory<C, P, O> &rs) = delete;
        calculator_factory &operator=(const calculator_factory<C, P, O> &rs) = delete;
        std::unordered_map<powerhouse::calculator_key, calculator_creator> _map;
        static std::once_flag only_one;
        static std::shared_ptr<calculator_factory<C, P, O>> _factory_instance;
    };
    template <typename C, typename P, typename O>
    std::shared_ptr<calculator_factory<C, P, O>> calculator_factory<C, P, O>::_factory_instance = nullptr;
    template <typename C, typename P, typename O>
    std::once_flag calculator_factory<C, P, O>::only_one;
}