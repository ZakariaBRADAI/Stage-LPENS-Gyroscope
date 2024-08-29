#define pinDIR 9
#define pinPUL 6

float omega = 2 * PI * 0.265;
float t, t0, dt, speed_RPM;
int N = 800;
bool mvt = false;

void setup() {
  Serial.begin(2000000);
  pinMode(pinPUL, OUTPUT);
  pinMode(pinDIR, OUTPUT);
}

void loop() {

  if (Serial.read() == '\n'){
    mvt = not(mvt);
    t0 = millis() / 1000.0;
  }

  if (mvt){
    t = millis() / 1000.0 - t0;

    speed_RPM =  17 * PI * omega * sin(omega * t) ;
    // Prefacteur à regler



    if  (speed_RPM != 0){
      dt = 1000 * 60 * 1000 / N / speed_RPM; 
      step(dt* 1000);
    }
  }
  
}


void step(float Dt){
  if (Dt>0) {
    digitalWrite(pinPUL, LOW);
    digitalWrite(pinDIR, HIGH);
    delayMicroseconds(Dt);

    digitalWrite(pinPUL, HIGH);
    digitalWrite(pinDIR, HIGH);
    delayMicroseconds(Dt);

    digitalWrite(pinPUL, HIGH);
    digitalWrite(pinDIR, LOW);
    delayMicroseconds(Dt);

    digitalWrite(pinPUL, LOW);
    digitalWrite(pinDIR, LOW);
    delayMicroseconds(Dt);
  }

  else if (Dt< 0){
    float dt2 = -Dt;
    digitalWrite(pinPUL, LOW);
    digitalWrite(pinDIR, LOW);
    delayMicroseconds(dt2);

    digitalWrite(pinPUL, HIGH);
    digitalWrite(pinDIR, LOW);
    delayMicroseconds(dt2);

    digitalWrite(pinPUL, HIGH);
    digitalWrite(pinDIR, HIGH);
    delayMicroseconds(dt2);

    digitalWrite(pinPUL, LOW);
    digitalWrite(pinDIR, HIGH);
    delayMicroseconds(dt2);
  }
  //return 0;
}
