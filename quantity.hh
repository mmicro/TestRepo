#include <chrono>
#include <functional>
#include <iostream>
#include <ratio>

#include <assert.h>
#include <string_view>
#include <type_traits>

// TODO: I will remove these lines at final version
#ifdef __clang__
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wc++98-compat"
    #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
    #pragma clang diagnostic ignored "-Wold-style-cast"
    #pragma clang diagnostic ignored "-Wunused-macros"
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wc++17-extensions"
#endif

namespace Optima::Control {

    namespace Prefix {
        // TODO: add hour, minute...
        using Atto = std::atto;
        using Femto = std::femto;
        using Pico = std::pico;
        using Nano = std::nano;
        using Micro = std::micro;
        using Milli = std::milli;
        using Centi = std::centi;
        using Deci = std::deci;
        using None = std::ratio<1, 1>;
        using Deca = std::deca;
        using Hecto = std::hecto;
        using Kilo = std::kilo;
        using Mega = std::mega;
        using Giga = std::giga;
        using Tera = std::tera;
        using Peta = std::peta;
        using Exa = std::exa;

        template <typename P1, typename P2>
        struct Multiply : std::ratio_multiply<P1, P2> {};

        template <typename P1, typename P2>
        struct Divide : std::ratio_divide<P1, P2> {};

        template <typename Ratio, std::size_t D>
        struct Pow {
            using type =
                typename Prefix::Multiply<Ratio,
                                          typename Pow<typename Ratio::type, D - 1>::type>::type;
        };

        template <typename Ratio>
        struct Pow<Ratio, 0> {
            using type = typename Ratio::type;
        };

        template <typename Ratio>
        struct Pow<Ratio, 1> {
            using type = typename Ratio::type;
        };

        template <typename P1, typename P2>
        constexpr inline bool IsEqual = std::ratio_equal_v<P1, P2>;
    }  // namespace Prefix

    namespace SI {

        using T = std::int32_t;
        static_assert(std::is_signed_v<T>, "SI underlying type must be signed");

#pragma region Primary_Dimension
        namespace details {
            // GCC has a problem now for implemnting class member specializatin in class! otherwise
            // I could make these lines less!
            template <T s, T m, T kg, T A, T K, T mol, T cd>
            constexpr inline T dim = 1;
            template <T s>
            constexpr inline T dim<s, 0, 0, 0, 0, 0, 0> = s < 0 ? -s : s;
            template <T m>
            constexpr inline T dim<0, m, 0, 0, 0, 0, 0> = m < 0 ? -m : m;
            template <T kg>
            constexpr inline T dim<0, 0, kg, 0, 0, 0, 0> = kg < 0 ? -kg : kg;
            template <T A>
            constexpr inline T dim<0, 0, 0, A, 0, 0, 0> = A < 0 ? -A : A;
            template <T K>
            constexpr inline T dim<0, 0, 0, 0, K, 0, 0> = K < 0 ? -K : K;
            template <T mol>
            constexpr inline T dim<0, 0, 0, 0, 0, mol, 0> = mol < 0 ? -mol : mol;
            template <T cd>
            constexpr inline T dim<0, 0, 0, 0, 0, 0, cd> = cd < 0 ? -cd : cd;
            // Scalar has no dimension
            template <>
            constexpr inline T dim<0, 0, 0, 0, 0, 0, 0> = 0;
        }  // namespace details
#pragma endregion

        /*
         * https://en.wikipedia.org/wiki/SI_base_unit
         */
        template <T s,    // time (second)
                  T m,    // length (meter)
                  T kg,   // mass (kilogram)
                  T A,    // current (ampere)
                  T K,    // thermodynamic temperature (kelvin)
                  T mol,  // amount of substance (mole)
                  T cd    // luminous intensity (candela)
                  >
        struct UnitSI {
            constexpr UnitSI operator+() const noexcept {
                return {};
            }
            constexpr UnitSI operator-() const noexcept {
                return {};
            }
            constexpr UnitSI operator+(UnitSI const&) const noexcept {
                return {};
            }
            constexpr UnitSI operator-(UnitSI const&) const noexcept {
                return {};
            }

            template <T os, T om, T okg, T oA, T oK, T oMol, T ocd>
            constexpr UnitSI<s + os, m + om, kg + okg, A + oA, K + oK, mol + oMol, cd + ocd>
            operator*(UnitSI<os, om, okg, oA, oK, oMol, ocd> const&) const noexcept {
                return {};
            }

            template <T os, T om, T okg, T oA, T oK, T oMol, T ocd>
            constexpr UnitSI<s - os, m - om, kg - okg, A - oA, K - oK, mol - oMol, cd - ocd>
            operator/(UnitSI<os, om, okg, oA, oK, oMol, ocd> const&) const noexcept {
                return {};
            }

            constexpr UnitSI& operator+=(UnitSI const&) noexcept {
                return *this;
            }

            constexpr UnitSI& operator-=(UnitSI const&) noexcept {
                return *this;
            }

            constexpr UnitSI& operator*=(UnitSI<0, 0, 0, 0, 0, 0, 0> const&) noexcept {
                return *this;
            }
            constexpr UnitSI& operator/=(UnitSI<0, 0, 0, 0, 0, 0, 0> const&) noexcept {
                return *this;
            }

            template <T os, T om, T okg, T oA, T oK, T oMol, T ocd>
            constexpr bool operator==(
                UnitSI<os, om, okg, oA, oK, oMol, ocd> const&) const noexcept {
                return s == os && m == om && kg == okg && A == oA && K == oK && mol == oMol &&
                       cd == ocd;
            }

            static constexpr auto dimension() noexcept {
                return details::dim<s, m, kg, A, K, mol, cd>;
            }
        };

#pragma region BASIC_SI
        using Scalar = UnitSI<0, 0, 0, 0, 0, 0, 0>;
        using Time = UnitSI<1, 0, 0, 0, 0, 0, 0>;
        using Length = UnitSI<0, 1, 0, 0, 0, 0, 0>;
        using Mass = UnitSI<0, 0, 1, 0, 0, 0, 0>;
        using Current = UnitSI<0, 0, 0, 1, 0, 0, 0>;
        using ThermoTemperature = UnitSI<0, 0, 0, 0, 1, 0, 0>;
        using SubstanceAmount = UnitSI<0, 0, 0, 0, 0, 1, 0>;
        using LuminousIntensity = UnitSI<0, 0, 0, 0, 0, 0, 1>;

        constexpr inline Scalar _scalar;
        constexpr inline Time _time;
        constexpr inline Length _length;
        constexpr inline Mass _mass;
        constexpr inline Current _current;
        constexpr inline ThermoTemperature _thremodynamic_temperature;
        constexpr inline SubstanceAmount _substance_amount;
        constexpr inline LuminousIntensity _luminous_intensity;
#pragma endregion

#pragma region DERIVED_QUANTITIES
        using Frequency = decltype(_scalar / _time);
        constexpr inline Frequency _frequency;

        using Velocity = decltype(_length / _time);
        constexpr inline Velocity _velocity;

        using Acceleration = decltype(_velocity / _time);
        constexpr inline Acceleration _acceleration;

        using Force = decltype(_mass * _acceleration);
        constexpr inline Force _force;

        using Power = decltype(_mass * _length * _length / (_time * _time * _time));
        constexpr inline Power _power;

        using Voltage = decltype(_power / _current);
        constexpr inline Voltage _voltage;

        using Resistance = decltype(_voltage / _current);
        constexpr inline Resistance _resistance;

        using Area = decltype(_length * _length);
        constexpr inline Area _area;

        using Volume = decltype(_area * _length);
        constexpr inline Volume _volume;
#pragma endregion

#pragma region VALIDATOR
        template <typename Q2>
        constexpr inline bool IsQuantityValid = false;
        template <>
        constexpr inline bool IsQuantityValid<SI::Scalar> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Time> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Length> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Mass> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Current> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::ThermoTemperature> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::SubstanceAmount> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::LuminousIntensity> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Frequency> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Velocity> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Acceleration> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Force> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Power> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Voltage> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Resistance> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Area> = true;
        template <>
        constexpr inline bool IsQuantityValid<SI::Volume> = true;
#pragma endregion

#pragma region UTILITY
        template <typename Q>
        constexpr inline auto quantityType = "scalar";
        template <>
        constexpr inline auto quantityType<SI::Time> = "time";
        template <>
        constexpr inline auto quantityType<SI::Length> = "length";
        template <>
        constexpr inline auto quantityType<SI::Mass> = "mass";
        template <>
        constexpr inline auto quantityType<SI::Current> = "current";
        template <>
        constexpr inline auto quantityType<SI::ThermoTemperature> = "thermodynamic_temperature";
        template <>
        constexpr inline auto quantityType<SI::SubstanceAmount> = "amount_of_substance";
        template <>
        constexpr inline auto quantityType<SI::LuminousIntensity> = "luminance";
        template <>
        constexpr inline auto quantityType<SI::Frequency> = "frequency";
        template <>
        constexpr inline auto quantityType<SI::Velocity> = "velocity";
        template <>
        constexpr inline auto quantityType<SI::Acceleration> = "acceleration";
        template <>
        constexpr inline auto quantityType<SI::Force> = "force";
        template <>
        constexpr inline auto quantityType<SI::Power> = "power";
        template <>
        constexpr inline auto quantityType<SI::Voltage> = "voltage";
        template <>
        constexpr inline auto quantityType<SI::Resistance> = "resistance";
        template <>
        constexpr inline auto quantityType<SI::Area> = "area";
        template <>
        constexpr inline auto quantityType<SI::Volume> = "volume";
#pragma endregion

    }  // namespace SI

    template <typename T, typename Q, typename P>
    class Quantity : private Q {
    public:
        using ValueType = T;
        using QType = Q;
        using Rep = P;

        static_assert(SI::IsQuantityValid<Q>, "Quantity type is not valid");

        Quantity() = default;
        explicit constexpr Quantity(T _value) noexcept : value{_value} {}
        explicit constexpr Quantity(T _value, Q) noexcept : value{_value} {}
        explicit constexpr Quantity(T _value, Q, P) noexcept : value{_value} {}

        Quantity(Quantity const&) = default;
        Quantity(Quantity&&) = default;

        Quantity& operator=(Quantity const&) = default;
        Quantity& operator=(Quantity&&) = default;

        constexpr ValueType Value() const noexcept {
            return value;
        }
        static constexpr QType Unit() noexcept {
            return {};
        }

        constexpr auto operator+() const noexcept {
            return Quantity{+value, +GetUnit(), Rep{}};
        }

        constexpr auto operator-() const noexcept {
            return Quantity{-value, -GetUnit(), Rep{}};
        }

        template <typename Q2, typename P2>
        constexpr auto operator+(Quantity<T, Q2, P2> const& other) const noexcept {
            using commonPrefix = std::common_type_t<Rep, P2>;
            return Quantity<T, Q2, commonPrefix>{convertToDst<commonPrefix, Rep>() +
                                                 other.template convertToDst<commonPrefix, P2>()};
        }

        template <typename Q2, typename P2>
        constexpr auto operator-(Quantity<T, Q2, P2> const& other) const noexcept {
            using commonPrefix = std::common_type_t<Rep, P2>;
            return Quantity<T, Q2, commonPrefix>{convertToDst<commonPrefix, Rep>() -
                                                 other.template convertToDst<commonPrefix, P2>()};
        }

        template <typename Q2, typename P2>
        constexpr auto operator*(Quantity<T, Q2, P2> const& other) const noexcept {
            auto resUnit = GetUnit() * other.GetUnit();
            using resQuantityType = decltype(resUnit);

            if constexpr (Q::dimension() > 1 || Q2::dimension() > 1 ||
                          resQuantityType::dimension() > 1) {
                using commonPrefix = std::common_type_t<P, P2>;
                return Quantity<T, resQuantityType, commonPrefix>{
                    convertToDst<commonPrefix, P>() *
                    other.template convertToDst<commonPrefix, P2>()};
            } else {
                using resRep = typename Prefix::Multiply<P, P2>::type;
                return Quantity<T, resQuantityType, resRep>{value * other.value};
            }
        }

        template <typename Q2, typename P2>
        constexpr auto operator/(Quantity<T, Q2, P2> const& other) const noexcept {
            auto resUnit = GetUnit() / other.GetUnit();
            using resQuantityType = decltype(resUnit);
            using commonPrefix = std::common_type_t<Rep, P2>;

            if constexpr (Q::dimension() > 1 || Q2::dimension() > 1 ||
                          resQuantityType::dimension() > 1) {
                return Quantity<T, resQuantityType, commonPrefix>{
                    convertToDst<commonPrefix, Rep>() /
                    other.template convertToDst<commonPrefix, P2>()};
            } else {
                using resRep = typename Prefix::Divide<P, P2>::type;
                return Quantity<T, resQuantityType, resRep>{value / other.value};
            }
        }

        // TODO: implement following operators
        // template <typename P2>
        // constexpr auto operator+=(Quantity<T, Q, P2, D> const& other) noexcept;
        // template <typename P2>
        // constexpr auto operator-=(Quantity<T, Q, P2, D> const& other) noexcept;
        // template <typename T2, typename Q2, typename P2, std::size_t D2>
        // constexpr auto operator*=(Quantity<T2, Q2, P2, D2> const& other) noexcept;
        // template <typename T2, typename Q2, typename P2, std::size_t D2>
        // constexpr auto operator/=(Quantity<T2, Q2, P2, D2> const& other) noexcept;

        constexpr auto getQuantityType() const noexcept {
            return std::string_view(SI::quantityType<Q>);
        }

        template <typename P2>
        constexpr auto cast() const noexcept {
            return Quantity<T, Q, P2>{this->convertToDst<P2, P>()};
        }

        template <typename Q2, typename P2>
        constexpr bool operator==(Quantity<T, Q2, P2> const& other) const noexcept {
            using comRep = std::common_type_t<Rep, P2>;
            return (std::is_same_v<Q, Q2>)&&(convertToDst<comRep, Rep>(),
                                             other.template convertToDst<comRep, P2>());
        }

        template <typename Q2, typename P2, typename Cmp = std::equal_to<>>
        constexpr bool isEqual(Quantity<T, Q2, P2> const& other, Cmp cmp = Cmp()) const noexcept {
            using comRep = std::common_type_t<Rep, P2>;
            return (std::is_same_v<Q, Q2>)&&(
                cmp(convertToDst<comRep, Rep>(), other.template convertToDst<comRep, P2>()));
        }

        static constexpr auto getPrefixRep() noexcept {
            return Rep{};
        }

        static constexpr auto getQuantityDim() noexcept {
            return Q::dimension();
        }

    private:
        template <typename T2, typename Q2, typename P2>
        friend class Quantity;

        ValueType value;

        constexpr QType& GetUnit() noexcept {
            return static_cast<QType&>(*this);
        }
        constexpr QType const& GetUnit() const noexcept {
            return static_cast<QType const&>(*this);
        }

        template <typename RepDst, typename RepSrc>
        constexpr T convertToDst() const noexcept {
            using Ratio = typename Prefix::Pow<typename std::ratio_divide<RepSrc, RepDst>::type,
                                               Q::dimension()>::type;
            return value * Ratio::num / Ratio::den;
        }

        constexpr T normalize() const noexcept {
            return value * Rep::num / Rep::den;
        }
    };

#pragma region QUANTITIES

#pragma region QUANTITY_SCALAR
    template <typename T = double>
    using Scalar = Quantity<T, SI::Scalar, Prefix::None>;
#pragma endregion

#pragma region QUANTITY_TIME
    template <typename T = double>
    using MegaSecond = Quantity<T, SI::Time, Prefix::Mega>;
    constexpr auto operator"" _mega_second(long double v) {
        return MegaSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mega_second(unsigned long long int v) {
        return MegaSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _MS(long double v) {
        return MegaSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _MS(unsigned long long int v) {
        return MegaSecond<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using KiloSecond = Quantity<T, SI::Time, Prefix::Kilo>;
    constexpr auto operator"" _kilo_second(long double v) {
        return KiloSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kilo_second(unsigned long long int v) {
        return KiloSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kS(long double v) {
        return KiloSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kS(unsigned long long int v) {
        return KiloSecond<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using Second = Quantity<T, SI::Time, Prefix::None>;
    constexpr auto operator"" _second(long double v) {
        return Second<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _second(unsigned long long int v) {
        return Second<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _S(long double v) {
        return Second<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _S(unsigned long long int v) {
        return Second<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MilliSecond = Quantity<T, SI::Time, Prefix::Milli>;
    constexpr auto operator"" _milli_second(long double v) {
        return MilliSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _milli_second(unsigned long long int v) {
        return MilliSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mS(long double v) {
        return MilliSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mS(unsigned long long int v) {
        return MilliSecond<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MicroSecond = Quantity<T, SI::Time, Prefix::Micro>;
    constexpr auto operator"" _micro_second(long double v) {
        return MicroSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _micro_second(unsigned long long int v) {
        return MicroSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uS(long double v) {
        return MicroSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uS(unsigned long long int v) {
        return MicroSecond<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using NanoSecond = Quantity<T, SI::Time, Prefix::Nano>;
    constexpr auto operator"" _nano_second(long double v) {
        return NanoSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nano_second(unsigned long long int v) {
        return NanoSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nS(long double v) {
        return NanoSecond<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nS(unsigned long long int v) {
        return NanoSecond<double>(static_cast<double>(v));
    }
#pragma endregion

#pragma region QUANTITY_CURRENT
    template <typename T = double>
    using MegaAmpere = Quantity<T, SI::Current, Prefix::Mega>;
    constexpr auto operator"" _mega_amp(long double v) {
        return MegaAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mega_amp(unsigned long long int v) {
        return MegaAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _MA(long double v) {
        return MegaAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _MA(unsigned long long int v) {
        return MegaAmpere<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using KiloAmpere = Quantity<T, SI::Current, Prefix::Kilo>;
    constexpr auto operator"" _kilo_amp(long double v) {
        return KiloAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kilo_amp(unsigned long long int v) {
        return KiloAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kA(long double v) {
        return KiloAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kA(unsigned long long int v) {
        return KiloAmpere<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using Ampere = Quantity<T, SI::Current, Prefix::None>;
    constexpr auto operator"" _amp(long double v) {
        return Ampere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _amp(unsigned long long int v) {
        return Ampere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _A(long double v) {
        return Ampere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _A(unsigned long long int v) {
        return Ampere<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MilliAmpere = Quantity<T, SI::Current, Prefix::Milli>;
    constexpr auto operator"" _milli_amp(long double v) {
        return MilliAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _milli_amp(unsigned long long int v) {
        return MilliAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mA(long double v) {
        return MilliAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mA(unsigned long long int v) {
        return MilliAmpere<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MicroAmpere = Quantity<T, SI::Current, Prefix::Micro>;
    constexpr auto operator"" _micro_amp(long double v) {
        return MicroAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _micro_amp(unsigned long long int v) {
        return MicroAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uA(long double v) {
        return MicroAmpere<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uA(unsigned long long int v) {
        return MicroAmpere<double>(static_cast<double>(v));
    }
#pragma endregion

#pragma region QUANTITY_VOLTAGE
    template <typename T = double>
    using KiloVolt = Quantity<T, SI::Voltage, Prefix::Kilo>;
    constexpr auto operator"" _kilo_volt(long double v) {
        return KiloVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kilo_volt(unsigned long long int v) {
        return KiloVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kV(long double v) {
        return KiloVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kV(unsigned long long int v) {
        return KiloVolt<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using Volt = Quantity<T, SI::Voltage, Prefix::None>;
    constexpr auto operator"" _volt(long double v) {
        return Volt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _volt(unsigned long long int v) {
        return Volt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _V(long double v) {
        return Volt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _V(unsigned long long int v) {
        return Volt<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MilliVolt = Quantity<T, SI::Voltage, Prefix::Milli>;
    constexpr auto operator"" _milli_volt(long double v) {
        return MilliVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _milli_volt(unsigned long long int v) {
        return MilliVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mV(long double v) {
        return MilliVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mV(unsigned long long int v) {
        return MilliVolt<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MicroVolt = Quantity<T, SI::Voltage, Prefix::Micro>;
    constexpr auto operator"" _micro_volt(long double v) {
        return MicroVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _micro_volt(unsigned long long int v) {
        return MicroVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uV(long double v) {
        return MicroVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uV(unsigned long long int v) {
        return MicroVolt<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using NanoVolt = Quantity<T, SI::Voltage, Prefix::Nano>;
    constexpr auto operator"" _nano_volt(long double v) {
        return NanoVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nano_volt(unsigned long long int v) {
        return NanoVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nV(long double v) {
        return NanoVolt<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _nV(unsigned long long int v) {
        return NanoVolt<double>(static_cast<double>(v));
    }
#pragma endregion

#pragma region QUANTITY_POWER
    template <typename T = double>
    using KiloWatt = Quantity<T, SI::Power, Prefix::Kilo>;

    template <typename T = double>
    using Watt = Quantity<T, SI::Power, Prefix::None>;

    template <typename T = double>
    using MilliWatt = Quantity<T, SI::Power, Prefix::Milli>;

    template <typename T = double>
    using MicroWatt = Quantity<T, SI::Power, Prefix::Micro>;
#pragma endregion

#pragma region QUANTITY_RESISTANCE
    template <typename T = double>
    using KiloOhm = Quantity<T, SI::Resistance, Prefix::Kilo>;
    constexpr auto operator"" _kilo_ohm(long double v) {
        return KiloOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kilo_ohm(unsigned long long int v) {
        return KiloOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kO(long double v) {
        return KiloOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _kO(unsigned long long int v) {
        return KiloOhm<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using Ohm = Quantity<T, SI::Resistance, Prefix::None>;
    constexpr auto operator"" _ohm(long double v) {
        return Ohm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _ohm(unsigned long long int v) {
        return Ohm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _O(long double v) {
        return Ohm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _O(unsigned long long int v) {
        return Ohm<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MilliOhm = Quantity<T, SI::Resistance, Prefix::Milli>;
    constexpr auto operator"" _milli_ohm(long double v) {
        return MilliOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _milli_ohm(unsigned long long int v) {
        return MilliOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mO(long double v) {
        return MilliOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _mO(unsigned long long int v) {
        return MilliOhm<double>(static_cast<double>(v));
    }

    template <typename T = double>
    using MicroOhm = Quantity<T, SI::Resistance, Prefix::Micro>;
    constexpr auto operator"" _micro_ohm(long double v) {
        return MicroOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _micro_ohm(unsigned long long int v) {
        return MicroOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uO(long double v) {
        return MicroOhm<double>(static_cast<double>(v));
    }
    constexpr auto operator"" _uO(unsigned long long int v) {
        return MicroOhm<double>(static_cast<double>(v));
    }
#pragma endregion

#pragma region QUANTITY_LENGTH
    template <typename T = double>
    using KiloMeter = Quantity<T, SI::Length, Prefix::Kilo>;

    template <typename T = double>
    using Meter = Quantity<T, SI::Length, Prefix::None>;

    template <typename T = double>
    using MilliMeter = Quantity<T, SI::Length, Prefix::Milli>;
#pragma endregion

#pragma region QUANTITY_MASS
    template <typename T = double>
    using KiloGram = Quantity<T, SI::Mass, Prefix::None>;

    template <typename T = double>
    using Gram = Quantity<T, SI::Mass, Prefix::Milli>;

    template <typename T = double>
    using MilliGram = Quantity<T, SI::Mass, Prefix::Micro>;
#pragma endregion

#pragma region QUANTITY_FORCE
    template <typename T = double>
    using KiloNewton = Quantity<T, SI::Force, Prefix::Kilo>;

    template <typename T = double>
    using Newton = Quantity<T, SI::Force, Prefix::None>;

    template <typename T = double>
    using MilliNewton = Quantity<T, SI::Force, Prefix::Milli>;
#pragma endregion

#pragma region QUANTITY_ACCELERATION
    template <typename T = double>
    using KiloMeterPerSecondSqr = Quantity<T, SI::Acceleration, Prefix::Kilo>;

    template <typename T = double>
    using MeterPerSecondSqr = Quantity<T, SI::Acceleration, Prefix::None>;

    template <typename T = double>
    using MilliMeterPerSecondSqr = Quantity<T, SI::Acceleration, Prefix::Milli>;
#pragma endregion

#pragma region QUANTITY_AREA
    template <typename T = double>
    using KiloMeterSqr = Quantity<T, SI::Area, Prefix::Kilo>;

    template <typename T = double>
    using MeterSqr = Quantity<T, SI::Area, Prefix::None>;

    template <typename T = double>
    using MilliMeterSqr = Quantity<T, SI::Area, Prefix::Milli>;
#pragma endregion

#pragma region QUANTITY_VOLUME
    template <typename T = double>
    using KiloMeterCub = Quantity<T, SI::Volume, Prefix::Kilo>;

    template <typename T = double>
    using MeterCub = Quantity<T, SI::Volume, Prefix::None>;

    template <typename T = double>
    using MilliMeterCub = Quantity<T, SI::Volume, Prefix::Milli>;
#pragma endregion

#pragma endregion

}  // namespace Optima::Control

#pragma region COMMON_TYPE
namespace std {
    using RatioType = decltype(micro::num);
    template <RatioType N0, RatioType D0, RatioType N1, RatioType D1>
    struct common_type<ratio<N0, D0>, ratio<N1, D1>> {
        using type = typename common_type_t<chrono::duration<int, ratio<N0, D0>>,
                                            chrono::duration<int, ratio<N1, D1>>>::period;
    };
}  // namespace std
#pragma endregion

#ifdef __clang__
    #pragma clang diagnostic pop
#endif
