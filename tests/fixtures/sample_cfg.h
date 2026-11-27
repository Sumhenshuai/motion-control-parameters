/* Synthetic header in the style of the Renesas RX24T motor library configuration.
 * Exercises comments, conditional compilation, numeric suffixes, hex literals and macro expressions. */
#ifndef SAMPLE_CFG_H
#define SAMPLE_CFG_H

#define MTR_SQRT_3                    (1.7320508f)
/*
#define MOTOR_SHUNT_TYPE_1_SHUNT
*/
#define MOTOR_SHUNT_TYPE_2_SHUNT      /* use two shunts */
#define MOTOR_COMMON_CFG_MTR_X        // selects the motor block below

#ifdef MOTOR_COMMON_CFG_MTR_X
#define MOTOR_CFG_POLE_PAIRS          (4)
#define MOTOR_CFG_MAGNETIC_FLUX       (0.0612f)
#define MOTOR_CFG_RESISTANCE          (0.85f)
#define MOTOR_CFG_D_INDUCTANCE        (0.0041f)
#define MOTOR_CFG_Q_INDUCTANCE        (0.0058f)
#define MOTOR_CFG_ROTOR_INERTIA       (1.2e-4f)
#define MOTOR_CFG_NOMINAL_CURRENT_RMS (6.0f)
#else
#define MOTOR_CFG_POLE_PAIRS          (99)
#endif

#define INVERTER_CFG_INPUT_V          (300.0f)
#if defined(MOTOR_SHUNT_TYPE_1_SHUNT)
#define MOTOR_MCU_CFG_CARRIER_FREQ    (10.0f)
#elif defined(MOTOR_SHUNT_TYPE_2_SHUNT)
#define MOTOR_MCU_CFG_CARRIER_FREQ    (8.0f)
#else
#define MOTOR_MCU_CFG_CARRIER_FREQ    (1.0f)
#endif
#define MOTOR_MCU_CFG_ADC_OFFSET      (0x7FF)
#define MOTOR_MCU_CFG_AD12BIT_DATA    (4095.0f)
#define MOTOR_COMMON_SPEED_LOOP_DIVIDER (4U)
#define CURRENT_CFG_PI_INTEGRAL_LIMIT_VQ (INVERTER_CFG_INPUT_V * 0.5f)
#define MOTOR_COMMON_PI_INTEGRAL_LIMIT_IQ (MOTOR_CFG_NOMINAL_CURRENT_RMS * MTR_SQRT_3)
#define CURRENT_CFG_OMEGA             (300.0f)
#define CURRENT_CFG_ZETA              (0.9f)
#define SPEED_CFG_OMEGA               (15.0f)
#define SPEED_CFG_ZETA                (1.0f)
#define SPEED_CFG_CTRL_PERIOD         (MOTOR_COMMON_SPEED_LOOP_DIVIDER / (MOTOR_MCU_CFG_CARRIER_FREQ * 1000))
#define CURRENT_CFG_FUNC_DUTY_SET     (R_Config_MTU3_MTU4_UpdDuty)
#endif /* SAMPLE_CFG_H */
