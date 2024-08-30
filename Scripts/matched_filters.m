triangle_pulse = [5 4 3 2 1] / sqrt(55); % triangluar pulse to be convolved with the bitstream.
bit_stream = randi([0, 1], 1, 10); %creates an array of 10 bits of zeroes and ones.
p_nrz_stream = 2*bit_stream - 1; %Turn the stream into polar nrz line code.
upsample_stream = upsample(p_nrz_stream, 5); %upsample so there's 5 bits every second.
signal_length = length(upsample_stream); % Calculate the length of the upsampled signal.
stream_time_ref = linspace(0, (signal_length - 1) * 0.2, signal_length);% Create an array of time reference for array of bits.

conv_length = length(triangle_pulse) + length(upsample_stream) - 1 ;
bit_stream_fourier = fft(upsample_stream, conv_length); %turn the bit polar nrz stream into frequncy domain equivalent
triangle_pulse_fourier = fft(triangle_pulse, conv_length); %turn the bit stream into frequncy domain equivalent
freq_domain_mult = bit_stream_fourier .* triangle_pulse_fourier; %multiply element wise
convolved_stream = ifft(freq_domain_mult,conv_length); % reverse fourier transfrom ( the convolved signal :D)
conv_stream_length = length(convolved_stream);
conv_time_ref = linspace(0, (conv_stream_length - 1) * 0.2, conv_stream_length);% Create an array of time reference for convolved stream.

figure;
stem(stream_time_ref, upsample_stream); % Plot the upsampled stream
xlabel('Time axis'); % X-axis label
ylabel('Amp axis'); % Y-axis label
title('Bit stream samples'); % Title

figure;
plot(triangle_pulse);
xlabel('Time axis'); % X-axis label
ylabel('Amp axis'); % Y-axis label
title('Triangle pulse shaper'); % Title

figure;
% 3.1 discrete plotting
stem(conv_time_ref, convolved_stream);
xlabel('Time axis'); % X-axis label
ylabel('Amp axis'); % Y-axis labels
title('Convolved output samples'); % Title

figure;
% 3.2 continous plotting
plot(convolved_stream);
xlabel('Time axis'); % X-axis label
ylabel('Amplitude axis'); % Y-axis labels
title('Convolved output'); % Title

matched_filter=fliplr(triangle_pulse);%%generating the matched filter
matchedfilter_output=conv(convolved_stream,matched_filter);
matchedfilter_time=linspace(0, (length(matchedfilter_output) - 1) * 0.2, length(matchedfilter_output));

rect=[1 1 1 1 1]/sqrt(5);
rect_output=conv(convolved_stream,rect);
rect_time=linspace(0, (length(rect_output) - 1) * 0.2, length(rect_output));

figure;
subplot(2,1,1);
plot(matchedfilter_time,matchedfilter_output);
hold on
stem(matchedfilter_time,matchedfilter_output,'blue');
xlabel('Time'); % X-axis label
ylabel('Amplitude'); % Y-axis label
title('output of The Matched Filter'); % Titlesubplot(2,1,2);
plot(rect_time,rect_output,'red');
hold on;
stem(rect_time,rect_output,'red');
xlabel('Time'); % X-axis label
ylabel('Amplitude'); % Y-axis label
title('output of The rectangular pulse'); % Title

for i=1:5:50
    sum=0;
    for j=0:1:4
        sum = sum + convolved_stream(i+j) * triangle_pulse(j+1);
        correlator_output(i+j) = sum;
    end
end

correlator_time=linspace(0, (length(correlator_output) - 1) * 0.2, length(correlator_output));

figure;
subplot(2,1,1);
plot(correlator_time,correlator_output,'red');
hold on;
stem(correlator_time,correlator_output,'red');
hold on
legend('Output of the correlator','output of the correlator (sampled)');
title("Correlator's Output"); % Title

subplot(2,1,2);
plot(matchedfilter_time,matchedfilter_output,'blue');
hold on
stem(matchedfilter_time,matchedfilter_output,'blue');
xlabel('Time'); % X-axis label
ylabel('Amplitude'); % Y-axis label
legend('Output of the matched filter','output of the matched filter (sampled)');
title("Matched Filter's Output"); % Title


bit_num = 10000;
data = randi([0 1], bit_num, 1); % Generates 10,000 random stream of bits
data_polar = (2 * data) - 1; % Maps the bit stream to 1 and -1 (Polar)
data_polar_Upsampled = upsample(data_polar, 5); % Upsample the stream

Eb = 1;

% Noise sould be added on the transmitted signal, that is convolved with
% the normalized pulse shaping function.

transmitted_signal = conv(data_polar_Upsampled, triangle_pulse);
BER_matchedFilter = zeros(1, 8);
BER_rect = zeros(1, 8);
BER_theoretical = zeros(1, 8);
miss_matchedFilter = 0;
miss_rectFilter = 0;

for SNR_db = -2 : 1 : 5
    SNR_linear = 10^(SNR_db / 10); % Linearize SNR to calculate N0
    N0 = Eb/SNR_linear;
    noise = sqrt(N0/2) * randn(length(transmitted_signal) , 1); % AWGN with variance = N0/2
    v = transmitted_signal + noise; % Insinuates the stream with the AWGN channel
    convolved_v_matchedFilter = conv(v, matched_filter);
    convolved_v_rect = conv(v, rect);

    % Calculating BER by comparing bits in the original signal with the
    % transmitted signal (that contains noise)

    for i = 1:1:10000
        if( ( data_polar(i) == 1 ) && ( convolved_v_matchedFilter(i*5) < 0 ) )
        
            miss_matchedFilter = miss_matchedFilter + 1;
        elseif( ( data_polar(i) == -1 ) && ( convolved_v_matchedFilter(i*5) > 0 ) )
        
            miss_matchedFilter = miss_matchedFilter + 1;
        end

        if( ( data_polar(i) == 1 ) && ( convolved_v_rect(i*5) < 0 ) )

            miss_rectFilter = miss_rectFilter + 1;

        elseif( ( data_polar(i) == -1 ) && ( convolved_v_rect(i*5) > 0 ) )

            miss_rectFilter = miss_rectFilter + 1;
        end
    end
 
    BER_matchedFilter(SNR_db +3) = miss_matchedFilter/ bit_num;
    BER_rect( SNR_db + 3) = miss_rectFilter / bit_num;

    BER_theoretical( SNR_db + 3) = 0.5 * erfc(sqrt(SNR_linear));

    miss_rectFilter = 0;
    miss_matchedFilter = 0;
end


figure;
SNR_db = -2:1:5;
semilogy(SNR_db, BER_matchedFilter, 'b-o'); % Plot BER_matchedFilter
hold on; % Hold the current plot

semilogy(SNR_db, BER_rect, 'r-*'); % Plot BER_rect

semilogy(SNR_db, BER_theoretical, 'g-'); % Plot BER_theoretical

xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate (BER) vs. SNR (semilog)');
legend('Matched Filter', 'Rect Filter', 'Theoretical'); % Add legend
grid on; % Enable grid
hold off; % Release the current plot

% Parameters
data_length = 100; % Length of transmitted data in bits

% Define the cases
cases = struct('R', [0, 0, 1, 1], 'delay', [2, 8, 2, 8]);
% Loop over cases

for i = 1:length(cases.R)
    R = cases.R(i);
    delay = cases.delay(i);
    
    % Generate filter coefficients using rcosdesign
    tx_filter = rcosdesign(R, delay, 5, 'sqrt');
    rx_filter = rcosdesign(R, delay, 5, 'sqrt');

    % Generate random data
    data = randi([0, 1], 1, data_length);

    % Convert 0's to -1's and 1's to +1's for PAM modulation
    PAM_SIGNAL = 2 * data - 1;

    % Edit PAM signal to have 5 samples to represent it
    Ts = 5; % Sampling period
    Edited_PAM_SIGNAL = upsample(PAM_SIGNAL, Ts);

    % Transmit through filter at point A
    transmitted_signal_A = filter( tx_filter, 1, Edited_PAM_SIGNAL);

    % Filter the transmitted signal at point B using the receive filter
    transmitted_signal_B = filter( rx_filter, 1, transmitted_signal_A);

    % Plot eye diagram for point A
    figure;
    eyediagram(transmitted_signal_A, 10 );
    title(sprintf('Eye Diagram for Point A (Case %d, R = %d, delay = %d)', i, R, delay));
    xlabel('Time');
    ylabel('Amplitude');

    % Plot eye diagram for point B
    figure;
    eyediagram(transmitted_signal_B, 10);
    title(sprintf('Eye Diagram for Point B (Case %d, R = %d, delay = %d)', i, R, delay));
    xlabel('Time');
    ylabel('Amplitude');
end



%% This section automates having to manually save each figure into a folder.

% Set the desired format (e.g., 'png', 'jpg', 'pdf', 'eps')
format = 'png';

% Get a list of all open figures
figureHandles = findall(0, 'type', 'figure');

% Iterate through each figure and save it
for i = 1:length(figureHandles)
    figure(figureHandles(i)); % Make the current figure active
    filename = sprintf('figure_%d.%s', i, format); % Generate a filename
    saveas(gcf, filename); % Save the figure
end