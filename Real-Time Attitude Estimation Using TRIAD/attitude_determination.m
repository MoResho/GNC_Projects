%% Clear Environment
clear;
clear mag;
clear imu;

%% Arduino and Sensor Initialization
a = arduino;
addr = '0x1E';                 % I2C address for QMC5883L
imu = mpu6050(a);              % Initialize MPU6050
mag = device(a, 'I2CAddress', addr);  % Initialize QMC5883L

% Configure HMC5883L Magnetometer
write(mag, [0x00 0x70]);  % Configuration Register A: 8-average, 15 Hz default, normal measurement
write(mag, [0x01 0x20]);  % Configuration Register B: Gain = 1.3 Ga
write(mag, [0x02 0x00]);  % Mode Register: Continuous measurement mode

pause(0.1);  % Small delay for sensor initialization

%% Reference Vectors in Navigation Frame
fn = [0.58, 0.30, 10.39];     % Reference acceleration vector
mn = [27.82, 4.64, -6.33];   % Reference magnetic vector
fn = fn / norm(fn);
mn = mn / norm(mn);

%% Main Loop
while true
    % Read sensor data
    fb = readAcceleration(imu);   % Accelerometer data (body frame)
    magVector = readMag(mag);     % Magnetometer data (body frame)
    mb = double([magVector(2), -magVector(1), magVector(3)]);  % Reorder & adjust signs

    % Normalize body-frame vectors
    W1 = fb / norm(fb);
    W2 = mb / norm(mb);

    % TRIAD algorithm - Body Frame Basis
    Ou1 = W1;
    Ou2 = cross(W1, W2) / norm(cross(W1, W2));
    Ou3 = cross(W1, Ou2) / norm(cross(W1, Ou2));
    Mou = [Ou1; Ou2; Ou3]';

    % TRIAD algorithm - Reference Frame Basis
    R1 = fn;
    R2 = cross(fn, mn) / norm(cross(fn, mn));
    R3 = cross(fn, R2) / norm(cross(fn, R2));
    Mr = [R1; R2; R3]';

    % Direction Cosine Matrix (DCM)
    Cbn = Mou * Mr';

    % Convert DCM to Quaternion and Plot
    q = ConvertAttitude(Cbn, 'DCM', 'EP');
    poseplot(quaternion(q'), MeshFileName = 'Satellite1.STL');
    drawnow;
end

%% Function: Read Magnetometer Data
function magData = readMag(mag)
    % Read 6 bytes from data output registers starting at 0x03
    write(mag, 0x03);  % Set pointer to data output X MSB
    rawData = read(mag, 6);

    % Extract signed 16-bit integers (Note: byte order is swapped)
    x = typecast(uint8([rawData(2), rawData(1)]), 'int16');
    z = typecast(uint8([rawData(4), rawData(3)]), 'int16');
    y = typecast(uint8([rawData(6), rawData(5)]), 'int16');

    % Convert raw data to microTesla (uT)
    x_uT = double(x) / 1100 * 100;
    y_uT = double(y) / 1100 * 100;
    z_uT = double(z) / 980  * 100;

    magData = [x_uT, y_uT, z_uT];
end
