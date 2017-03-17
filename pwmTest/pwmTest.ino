// PWM Read/Write Module

/* Kalman Filter - What's missing:
 * Motor constants
 * Initial state
 * Initial error covariance
 * Process noise covariance
 * Measurement noise covariance
 *
 * 
 * Resources:
 * https://www.pololu.com/product/1444/faqs
 * 
 * Matlab Commands
 * delete(instrfindall) = clear all com ports
 * instrfind = print all com ports
 * 
 */

#define MATLAB_OUTPUT

#define HallPin1 2
#define HallPin2 3
#define PwmPin1 13
#define PwmPin2 5

#define SAMPLE_COUNT_MAX_HALL 20

#define FRICTION   0.11                  // b
#define MOMENT     0.02                  // J
#define TORQUE     0.01           // K = stall torque / stall current
#define INDUCTANCE .5                  // L
#define RESISTANCE 2.7                // R (measured across the motor with ohmmeter)
                                      // Voltage Constant = 0.92765 = 4.08v / 0.7Hz (in rad/s)
                                      // Torque = Voltage Constant / 0.011827 


volatile unsigned int dataReadyHall;
char incomingByte = 0;
char incomingArray[8];
long previousTime = 0;
long deltaTimeHall = 0;
double avgPeriodHall =0;
double avgFrequency;
unsigned char sampleCountHall = 0;
volatile long encoderCount =0;
long printTimeStamp;


const byte numChars = 32;
char receivedChars[numChars]; // an array to store the received data
signed int receivedInt;

boolean newData = false;

// system matrix
double a00 = (- FRICTION / MOMENT);
double a01 = (TORQUE / MOMENT);
double a10 = (-TORQUE / INDUCTANCE);
double a11 = (- RESISTANCE / INDUCTANCE);

// state
double x0 = 0.2;
double x1 = 0.3;

// state prediction
double xP0 = 0.2;
double xP1 = 0.3;

// measurement
double z0 = 0;
double z1 = 0;

// process noise covariance
double q00 = 1;
double q01 = 0;
double q10 = 0;
double q11 = q00;


// error covariance
double p00 = 10*q00;
double p01 = 0;
double p10 = 0;
double p11 = 10*q11;

// error covariance prediction
double pP00 = 10*q00;
double pP01 = 0;
double pP10 = 0;
double pP11 = 10*q11;

// measurement equation matrix
double h00 = 1;
double h01 = 0;
//double h10 = 0;
//double h11 = 0;



// measurement noise covariance
double r00 = 0.02;//0.01672098793;
double r01 = 0;
double r10 = 0;
double r11 = 1;  // returns NaN if this is 0 because of detHPH div0 error

double voltage = 0;

void kalman_TimeUpdate (void)
{
  // Predict the next state based off of the previous state.
  xP0 = a00 * x0 + a01 * x1;
  xP1 = a10 * x0 + a11 * x1 + ( 1 / INDUCTANCE) * voltage;

  printData("x1: ", x1, 2);
//----------------------------------------------
  // Predict the next error covariance based off of the previous error covariance

  // A * P
  double ap00, ap01, ap10, ap11;
  
  ap00 = a00 * p00 + a01 * p10;
  ap01 = a00 * p01 + a01 * p11;
  ap10 = a10 * p00 + a11 * p10;
  ap11 = a10 * p01 + a11 * p11;

  // (A * P) * A^T
  pP00 = ap00 * a00 + ap01 * a01;
  pP01 = ap00 * a10 + ap01 * a11;
  pP10 = ap10 * a00 + ap11 * a01;
  pP11 = ap10 * a10 + ap11 * a11;

  //printData("pP00: ", pP00, 2);

  // [(A*P)*A^T] + Q
  pP00 += q00;
  pP01 += q01;
  pP10 += q10;
  pP11 += q11;
}

void kalman_MeasurementUpdateTest (void)
{
  // Compute the Kalman Gain

  // H * P
  double hp00, hp01, hp10, hp11;
  
  hp00 = h00 * pP00 + h01 * pP10;
  hp01 = h00 * pP01 + h01 * pP11;

  // (H*P) * H^T
  double hph00, hph01, hph10, hph11;
  
  hph00 = hp00 * h00 + hp01 * h01;

  // [(H*P)*H^T] + R
  hph00 += r00;

  // [[(H*P)*H^T]+R]^-1
  double hphi00, hphi01, hphi10, hphi11;

  hphi00 = hph00;

  // P * H^T
  double ph00, ph01, ph10, ph11;
  
  ph00 = pP00 * h00 + pP01 * h01;
  ph10 = pP10 * h00 + pP11 * h01;
  
  // (P*H^T) * [[(H*P)*H^T]+R]^-1
  double k00, k01, k10, k11;
  
  k00 = ph00 / hphi00;
  k10 = ph10 / hphi00;


//----------------------------------------------

  // Update the state estimate with measurement

  // (H * x)
  double hx0, hx1;
  hx0 = h00 * xP0 + h01 * xP1;
  //hx1 = h10 * xP0 + h11 * xP1;

  // z - (H*x)
  double zhx0, zhx1;
  zhx0 = z0 - hx0;
  //zhx1 = z1 - hx1;

  // K * [z-(H*x)]
  double kzhx0, kzhx1;
  kzhx0 = k00 * zhx0;
  kzhx1 = k10 * zhx0;

  // x + [K*[z-(H*x)]]
  x0 = xP0 + kzhx0;
  x1 = xP1 + kzhx1;

  printData("zhx0: ", zhx0, 2);
//----------------------------------------------

  // Update the error covariance

  // K * H
  double kh00, kh01, kh10, kh11;
  
  kh00 = k00 * h00;
  kh01 = k00 * h01;
  kh10 = k10 * h00;
  kh11 = k10 * h01;

  // I - (K*H)
  double ikh00, ikh01, ikh10, ikh11;
  
  ikh00 = 1 - kh00;
  ikh01 = 0 - kh01;
  ikh10 = 0 - kh10;
  ikh11 = 1 - kh11;

  // [I-(K*H)] * P
  p00 = ikh00 * pP00 + ikh01 * pP10;
  p01 = ikh00 * pP01 + ikh01 * pP11;
  p10 = ikh10 * pP00 + ikh11 * pP10;
  p11 = ikh10 * pP01 + ikh11 * pP11; 

}

/*
void kalman_MeasurementUpdate (void)
{
  // Compute the Kalman Gain

  // H * P
  double hp00, hp01, hp10, hp11;
  
  hp00 = h00 * pP00 + h01 * pP01;
  hp01 = h00 * pP10 + h01 * pP11;
  hp10 = h10 * pP00 + h11 * pP01;
  hp11 = h10 * pP10 + h11 * pP11;

  printData("pP00: ", pP00, 2);

  // (H*P) * H^T
  double hph00, hph01, hph10, hph11;
  
  hph00 = hp00 * h00 + hp01 * h10;
  hph01 = hp00 * h01 + hp01 * h11;
  hph10 = hp10 * h00 + hp11 * h10;
  hph11 = hp10 * h01 + hp11 * h11;

  // [(H*P)*H^T] + R
  hph00 += r00;
  hph01 += r01;
  hph10 += r10;
  hph11 += r11;

//  printData("hph00: ", hph00, 2);
//  printData("hph01: ", hph01, 2);
//  printData("hph10: ", hph10, 2);
//  printData("hph11: ", hph11, 2);
    
  // [[(H*P)*H^T]+R]^-1
  double detHPH = hph00 * hph11 - hph01 * hph10;
  double hphi00, hphi01, hphi10, hphi11;

  printData("detHPH: ", detHPH, 2);

  if(detHPH == 0)
  {
    detHPH = 1;
  }
  
  hphi00 = hph11 / detHPH;
  hphi01 = -hph01 / detHPH;
  hphi10 = -hph10 / detHPH;
  hphi11 = hph00 / detHPH;

  
  
  

  // P * H^T
  double ph00, ph01, ph10, ph11;
  
  ph00 = pP00 * h00 + pP01 * h10;
  ph01 = pP00 * h01 + pP01 * h11;
  ph10 = pP10 * h00 + pP11 * h10;
  ph11 = pP10 * h01 + pP11 * h11;

  // (P*H^T) * [[(H*P)*H^T]+R]^-1
  double k00, k01, k10, k11;
  
  k00 = ph00 * hphi00 + ph01 * hphi01;
  k01 = ph00 * hphi10 + ph01 * hphi11;
  k10 = ph10 * hphi00 + ph11 * hphi01;
  k11 = ph10 * hphi10 + ph11 * hphi11;


//----------------------------------------------

  // Update the state estimate with measurement

  // (H * x)
  double hx0, hx1;
  hx0 = h00 * xP0 + h01 * xP1;
  hx1 = h10 * xP0 + h11 * xP1;

  // z - (H*x)
  double zhx0, zhx1;
  zhx0 = z0 - hx0;
  zhx1 = z1 - hx1;

  // K * [z-(H*x)]
  double kzhx0, kzhx1;
  kzhx0 = k00 * zhx0 + k01 * zhx1;
  kzhx1 = k10 * zhx0 + k11 * zhx1;

  // x + [K*[z-(H*x)]]
  x0 = xP0 + kzhx0;
  x1 = xP1 + kzhx1;
//----------------------------------------------

  // Update the error covariance

  // K * H
  double kh00, kh01, kh10, kh11;
  
  kh00 = k00 * h00 + k01 * h01;
  kh01 = k00 * h10 + k01 * h11;
  kh10 = k10 * h00 + k11 * h01;
  kh11 = k10 * h10 + k11 * h11;

  // I - (K*H)
  double ikh00, ikh01, ikh10, ikh11;
  
  ikh00 = 1 - kh00;
  ikh01 = 0 - kh01;
  ikh10 = 0 - kh10;
  ikh11 = 1 - kh11;

  // [I-(K*H)] * P
  p00 = ikh00 * pP00 + ikh01 * pP01;
  p01 = ikh00 * pP10 + ikh01 * pP11;
  p10 = ikh10 * pP00 + ikh11 * pP01;
  p11 = ikh10 * pP10 + ikh11 * pP11; 
}

*/

void printData (String Name, double Data, unsigned char DecimalPlaces)
{
#ifndef MATLAB_OUTPUT
  Serial.print(Name);
  Serial.print(Data, DecimalPlaces);
  Serial.print("\t");
 #endif
}

void setup() 
{   
  Serial.begin(115200);   // Start serial connection

  pinMode(PwmPin1, OUTPUT);
  pinMode(HallPin1, INPUT);
  setPwmDutyCycle(PwmPin1, 60);
  pinMode(HallPin2, INPUT);
  attachInterrupt(digitalPinToInterrupt(HallPin1), isr_readPwm, RISING);

  // Declare the PWM pins
  
  previousTime = micros();
}

void loop() 
{
  recvWithEndMarker();
  showNewData();

  if(dataReadyHall)
  { 
    readHallData();
        
    z0 = avgFrequency; //* 6.2831853; // hz to rad/s

    //if(millis()>5000 && receivedInt != 0)
    {
      kalman_TimeUpdate();
      kalman_MeasurementUpdateTest();
    

      printData("Frequency: ", avgPeriodHall, 4);
      printData("RecievedInt: ", receivedInt, 0);
      printData("x0: ", x0, 2);
      printData("xP0: ", xP0, 2);

      if(micros() - printTimeStamp > 5)
      {
        printTimeStamp = micros();
        Serial.print(x0,6);
        //Serial.print("/t");
        Serial.println();
        Serial.print(avgPeriodHall,6);
    
        Serial.println();
      }
    }
    dataReadyHall = 0;
  }
  
  if(newData)
  {
  #ifndef MATLAB_OUTPUT
    Serial.println("Recieved new duty cycle...");
  #else
    receivedInt = 60;
  #endif
  
    voltage = 2 + 3 * ((double)receivedInt/100);
    // 60% = 3.85V
    // 30% = 2.75V

    
    setPwmDutyCycle(PwmPin1, receivedInt);
    newData = false;
  }

  
}





void readHallData(void)
{
  deltaTimeHall = micros() - previousTime;
  previousTime = micros();
  
  avgPeriodHall = (double) encoderCount / deltaTimeHall;   // encoder ticks / microsecond

  encoderCount = 0;
  
  avgPeriodHall *= 1000000;                                     // encoder ticks / second
  avgPeriodHall /= 16;                                    // revolutions / second (or Hertz)
  avgPeriodHall /= 50;

  if(sampleCountHall < SAMPLE_COUNT_MAX_HALL)
  {
    sampleCountHall++;
  }
  
  avgFrequency += ((avgPeriodHall/(1)) - avgFrequency) / sampleCountHall;

}


// @param 'pin'       : the PWM pin we are setting
// @param 'dutyCycle' : the duty cycle we wish to set it at, integer from -100% to +100%
void setPwmDutyCycle(const unsigned int pin, signed int dutyCycle)
{

  float dutyCycleTemp = dutyCycle;

  // convert to decimal
  dutyCycleTemp /= 100;
  dutyCycleTemp += 1;
  dutyCycleTemp *= 127.5;

  // round to the nearest integer
  if (dutyCycleTemp > 127)
  {
    dutyCycle = (signed int)ceil(dutyCycleTemp);
  }
  else
  {
    dutyCycle = (signed int)floor(dutyCycleTemp);
  }

  //  set the pin accordingly
  analogWrite(pin, dutyCycle);
}

void isr_readPwm (void)
{
  encoderCount++;
  
  if(dataReadyHall == 0)
  {
    dataReadyHall = 1;
  }
}

void recvWithEndMarker() 
{
 static byte ndx = 0;
 char endMarker = '\n';
 char rc;
 
 while (Serial.available() > 0 && newData == false) 
 {
   rc = Serial.read();
  
   if (rc != endMarker) 
   {
     receivedChars[ndx] = rc;
     ndx++;
   
     if (ndx >= numChars) 
     {
      ndx = numChars - 1;
     }
   }
   else 
   {
     receivedChars[ndx] = '\0'; // terminate the string
     ndx = 0;
     newData = true;

     receivedInt = atoi(receivedChars);
   }
 }
}

void showNewData() 
{
#ifndef MATLAB_OUTPUT
 if (newData == true) 
 {
   Serial.print("Serial Input:  ");
   Serial.println(receivedInt);
   //newData = false;
 }
#endif
}

