/*
 * This file is part of libsidplayfp, a SID player engine.
 *
 * Copyright 2011-2019 Leandro Nini <drfiemost@users.sourceforge.net>
 * Copyright 2007-2010 Antti Lankila
 * Copyright 2010 Dag Lem
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include "FilterModelConfig.h"

#ifdef __LIBRETRO__
#ifdef __PS3__
#include "PS3_include.h"
#endif
#else
#include <cmath>
#include <cassert>
#endif

#include "Integrator.h"
#include "OpAmp.h"

namespace reSIDfp
{

#ifndef HAVE_CXX11
inline double log1p(double x)
{
    return log(1. + x) - (((1. + x) - 1.) - x) / (1. + x);
}
#endif

const unsigned int OPAMP_SIZE = 33;

const Spline::Point opamp_voltage[OPAMP_SIZE] =
{
  {  0.81, 10.31 },
  {  2.40, 10.31 },
  {  2.60, 10.30 },
  {  2.70, 10.29 },
  {  2.80, 10.26 },
  {  2.90, 10.17 },
  {  3.00, 10.04 },
  {  3.10,  9.83 },
  {  3.20,  9.58 },
  {  3.30,  9.32 },
  {  3.50,  8.69 },
  {  3.70,  8.00 },
  {  4.00,  6.89 },
  {  4.40,  5.21 },
  {  4.54,  4.54 },
  {  4.60,  4.19 },
  {  4.80,  3.00 },
  {  4.90,  2.30 },
  {  4.95,  2.03 },
  {  5.00,  1.88 },
  {  5.05,  1.77 },
  {  5.10,  1.69 },
  {  5.20,  1.58 },
  {  5.40,  1.44 },
  {  5.60,  1.33 },
  {  5.80,  1.26 },
  {  6.00,  1.21 },
  {  6.40,  1.12 },
  {  7.00,  1.02 },
  {  7.50,  0.97 },
  {  8.50,  0.89 },
  { 10.00,  0.81 },
  { 10.31,  0.81 },
};

std::unique_ptr<FilterModelConfig> FilterModelConfig::instance(nullptr);

FilterModelConfig* FilterModelConfig::getInstance()
{
    if (!instance.get())
    {
        instance.reset(new FilterModelConfig());
    }
    return instance.get();
}

void FilterModelConfig::initTables()
{
    dac.kinkedDac(MOS6581);

    Spline::Point scaled_voltage[OPAMP_SIZE];

    for (unsigned int i = 0; i < OPAMP_SIZE; i++)
    {
        scaled_voltage[i].x = N16 * (opamp_voltage[i].x - opamp_voltage[i].y + denorm) / 2.;
        scaled_voltage[i].y = N16 * (opamp_voltage[i].x - vmin);
    }

    Spline s(scaled_voltage, OPAMP_SIZE);

    for (int x = 0; x < (1 << 16); x++)
    {
        const Spline::Point out = s.evaluate(x);
        double tmp = out.x;
        if (tmp < 0.) tmp = 0.;
        assert(tmp < 65535.5);
        opamp_rev[x] = static_cast<unsigned short>(tmp + 0.5);
    }

    OpAmp opampModel(opamp_voltage, OPAMP_SIZE, kVddt);

    for (int i = 0; i < 5; i++)
    {
        const int idiv = 2 + i;
        const int size = idiv << 16;
        const double n = idiv;
        opampModel.reset();
        summer[i] = new unsigned short[size];

        for (int vi = 0; vi < size; vi++)
        {
            const double vin = vmin + vi / N16 / idiv;
            const double tmp = (opampModel.solve(n, vin) - vmin) * N16;
            assert(tmp > -0.5 && tmp < 65535.5);
            summer[i][vi] = static_cast<unsigned short>(tmp + 0.5);
        }
    }

    for (int i = 0; i < 8; i++)
    {
        const int idiv = (i == 0) ? 1 : i;
        const int size = (i == 0) ? 1 : i << 16;
        const double n = i * 8.0 / 6.0;
        opampModel.reset();
        mixer[i] = new unsigned short[size];

        for (int vi = 0; vi < size; vi++)
        {
            const double vin = vmin + vi / N16 / idiv;
            const double tmp = (opampModel.solve(n, vin) - vmin) * N16;
            assert(tmp > -0.5 && tmp < 65535.5);
            mixer[i][vi] = static_cast<unsigned short>(tmp + 0.5);
        }
    }

    for (int n8 = 0; n8 < 16; n8++)
    {
        const int size = 1 << 16;
        const double n = n8 / 8.0;
        opampModel.reset();
        gain[n8] = new unsigned short[size];

        for (int vi = 0; vi < size; vi++)
        {
            const double vin = vmin + vi / N16;
            const double tmp = (opampModel.solve(n, vin) - vmin) * N16;
            assert(tmp > -0.5 && tmp < 65535.5);
            gain[n8][vi] = static_cast<unsigned short>(tmp + 0.5);
        }
    }

    const double nkVddt = N16 * kVddt;
    const double nVmin = N16 * vmin;

    for (unsigned int i = 0; i < (1 << 16); i++)
    {
        const double Vg = nkVddt - sqrt((double)(i << 16));
        const double tmp = k * Vg - nVmin;
        assert(tmp > -0.5 && tmp < 65535.5);
        vcr_kVg[i] = static_cast<unsigned short>(tmp + 0.5);
    }

    const double kVt = k * Vth;
    const double Is = 2. * uCox * Ut * Ut / k * WL_vcr;
    const double N15 = norm * ((1 << 15) - 1);
    const double n_Is = N15 * 1.0e-6 / C * Is;

    for (int kVg_Vx = 0; kVg_Vx < (1 << 16); kVg_Vx++)
    {
        const double log_term = log1p(exp((kVg_Vx / N16 - kVt) / (2. * Ut)));
        const double tmp = n_Is * log_term * log_term;
        assert(tmp > -0.5 && tmp < 65535.5);
        vcr_n_Ids_term[kVg_Vx] = static_cast<unsigned short>(tmp + 0.5);
    }
}

FilterModelConfig::FilterModelConfig() :
    voice_voltage_range(1.5),
    voice_DC_voltage(5.0),
    C(470e-12),
    Vdd(12.18),
    Vth(1.31),
    Ut(26.0e-3),
    k(1.0),
    uCox(20e-6),
    WL_vcr(9.0 / 1.0),
    WL_snake(1.0 / 115.0),
    kVddt(k * (Vdd - Vth)),
    dac_zero(6.65),
    dac_scale(2.63),
    vmin(opamp_voltage[0].x),
    vmax(kVddt < opamp_voltage[0].y ? opamp_voltage[0].y : kVddt),
    denorm(vmax - vmin),
    norm(1.0 / denorm),
    N16(norm * ((1 << 16) - 1)),
    dac(DAC_BITS)
{
    initTables();
}

FilterModelConfig::~FilterModelConfig()
{
    for (int i = 0; i < 5; i++)
    {
        delete [] summer[i];
    }
    for (int i = 0; i < 8; i++)
    {
        delete [] mixer[i];
    }
    for (int i = 0; i < 16; i++)
    {
        delete [] gain[i];
    }
}

unsigned short* FilterModelConfig::getDAC(double adjustment) const
{
    const double dac_zero = getDacZero(adjustment);
    unsigned short* f0_dac = new unsigned short[1 << DAC_BITS];

    for (unsigned int i = 0; i < (1 << DAC_BITS); i++)
    {
        const double fcd = dac.getOutput(i);
        const double tmp = N16 * (dac_zero + fcd * dac_scale / (1 << DAC_BITS) - vmin);
        assert(tmp > -0.5 && tmp < 65535.5);
        f0_dac[i] = static_cast<unsigned short>(tmp + 0.5);
    }

    return f0_dac;
}

std::unique_ptr<Integrator> FilterModelConfig::buildIntegrator()
{
    double tmp = N16 * (kVddt - vmin);
    assert(tmp > -0.5 && tmp < 65535.5);
    const unsigned short nkVddt = static_cast<unsigned short>(tmp + 0.5);

    tmp = denorm * (1 << 13) * (uCox / (2. * k) * WL_snake * 1.0e-6 / C);
    assert(tmp > -0.5 && tmp < 65535.5);
    const unsigned short n_snake = static_cast<unsigned short>(tmp + 0.5);

    return std::unique_ptr<Integrator>(new Integrator(vcr_kVg, vcr_n_Ids_term, opamp_rev, nkVddt, n_snake));
}

} // namespace reSIDfp