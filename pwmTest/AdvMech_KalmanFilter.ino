/* Kalman Filter for AdvMech Lab 3
 * https://github.com/Krelkal/AdvMech-Lab3
 *  
 *  
 *  
 */

//---------------------------------------------------------------------------
//-- Defines
//---------------------------------------------------------------------------

#define MATLAB_OUTPUT
#define MATLAB_STATIC_DUTY_CYCLE 60

#define HallPin1 2
#define HallPin2 3
#define PwmPin1 13
#define PwmPin2 5

#define SAMPLE_COUNT_MAX_HALL 20

#define FRICTION   0.11   // b
#define MOMENT     0.02   // J
#define TORQUE     0.01   // K 
#define INDUCTANCE .5     // L
#define RESISTANCE 2.7    // R

//---------------------------------------------------------------------------
//-- Variables
//---------------------------------------------------------------------------

// isr
volatile unsigned int dataReadyHall;
volatile long encoderCount =0;

// hall effect sensor
long previousTime = 0;
long deltaTimeHall = 0;
double frequency = 0;
double avgFrequency;
unsigned char sampleCountHall = 0;

// serial communication
const byte numChars = 32;
char receivedChars[numChars]; // an array to store the received data
signed int receivedInt;
boolean newSerialData = false;

// misc
long printTimeStamp;
double voltage = 0;

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

// measurement noise covariance
double r00 = 0.02;//0.01672098793;

//---------------------------------------------------------------------------
//-- Functions
//---------------------------------------------------------------------------

void setup() 
{
  // Start serial connection 
  Serial.begin(115200);   

  // Init pins
  pinMode(PwmPin1, OUTPUT);
  pinMode(HallPin1, INPUT);
  pinMode(HallPin2, INPUT);
  attachInterrupt(digitalPinToInterrupt(HallPin1), isr_readPwm, RISING);
  
#ifdef MATLAP_OUTPUT
  // Set initial speed for MATLAB output
  setPwmDutyCycle(PwmPin1, MATLAB_STATIC_DUTY_CYCLE);
#endif

  // First time stamp
  previousTime = micros();
}

void loop() 
{
  // Watch for incoming command
  heartbeat_readSerial();
  
  #ifndef MATLAB_OUTPUT
    heartbeat_echoNewData();
  #endif

  if(dataReadyHall) // isr flag
  { 
    readHall();

    // uh ohhh... should this be frequency or avgFrequency?!
    z0 = avgFrequency; //* 6.2831853; // hz to rad/s

#ifndef MATLAB_OUTPUT
    // wait five seconds and for a non-zero input before starting the motor (safety yo)
    if(millis()>5000 && receivedInt != 0)
#endif
    {
      // run Kalman filter
      kalman_TimeUpdate();
      kalman_MeasurementUpdate();
    
      // print useful data and stuff
      printData("Frequency: ", frequency, 4);
      printData("RecievedInt: ", receivedInt, 0);
      printData("x0: ", x0, 2);
      printData("xP0: ", xP0, 2);

      if(micros() - printTimeStamp > 5) // don't print every itteration
      {
        printTimeStamp = micros();
        Serial.print(x0,6);
        Serial.println();
        Serial.print(frequency, 6);
        Serial.println();
      }
    }
    dataReadyHall = 0;
  }

  // update the motor speed
  if(newSerialData)
  {
    newSerialData = false;
    
  #ifndef MATLAB_OUTPUT
    Serial.println("Recieved new duty cycle...");
  #else
    
    receivedInt = MATLAB_STATIC_DUTY_CYCLE;
  #endif

    // 'estimate' the voltage based on the duty cycle
    voltage = 2 + 3 * ((double)receivedInt/100);

    setPwmDutyCycle(PwmPin1, receivedInt);
  }
}

void kalman_TimeUpdate (void)
{
  // Predict the next state based off of the previous state.
  xP0 = a00 * x0 + a01 * x1;
  xP1 = a10 * x0 + a11 * x1 + ( 1 / INDUCTANCE) * voltage;

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

  // [(A*P)*A^T] + Q
  pP00 += q00;
  pP01 += q01;
  pP10 += q10;
  pP11 += q11;
}

void kalman_MeasurementUpdate (void)
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

// @param 'Name'          : a string that will be printed before the data
// @param 'Data'          : the data to be printed
// @param 'DecimalPlaces' : how many decimal places will be displayed
void printData (String Name, double Data, unsigned char DecimalPlaces)
{
#ifndef MATLAB_OUTPUT
  Serial.print(Name);
  Serial.print(Data, DecimalPlaces);
  Serial.print("\t");
 #endif
}

void readHall(void)
{
  // Find the delta time since last read
  deltaTimeHall = micros() - previousTime;
  previousTime = micros();

  
  frequency = (double) encoderCount / deltaTimeHall;   // encoder ticks / microsecond

  encoderCount = 0; // reset encoder counter
  
  frequency *= 1000000;  // encoder ticks / second
  frequency /= 16;       // revolutions / second of encoder
  frequency /= 50;       // revolutions / second of motor shaft

  // increment counter
  if(sampleCountHall < SAMPLE_COUNT_MAX_HALL)
  {
    sampleCountHall++;
  }

  // low pass filter
  avgFrequency += (frequency- avgFrequency) / sampleCountHall;

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

void heartbeat_readSerial(void) 
{
 static byte ndx = 0;
 char endMarker = '\n';
 char rc;
 
 while (Serial.available() > 0 && newSerialData == false) 
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
     newSerialData = true;
     receivedInt = atoi(receivedChars);
   }
 }
}

void heartbeat_echoNewData(void) 
{
 if (newSerialData == true) 
 {
   Serial.print("Serial Input:  ");
   Serial.println(receivedInt);
 }
}

void isr_readPwm (void)
{
  // increment counter
  encoderCount++;

  // set flag
  if(dataReadyHall == 0)
  {
    dataReadyHall = 1;
  }
}
