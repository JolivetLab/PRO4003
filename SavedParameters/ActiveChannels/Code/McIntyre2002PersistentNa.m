function activeChannel = McIntyre2002PersistentNa()

activeChannel.channames =                           'Persistent Na+ (McIntyre 2002)';
activeChannel.cond.value.ref =                      0.05;
activeChannel.cond.value.vec =                      [];
activeChannel.cond.units =                          {2, 'mS', 'mm', [1, -2]};
activeChannel.erev.value =                          60;
activeChannel.erev.units =                          {1, 'mV', 1};
activeChannel.gates.temp =                          36;
activeChannel.gates.number =                        1;
activeChannel.gates.label =                         {'p'};
activeChannel.gates.numbereach =                    3;
activeChannel.gates.alpha.q10 =                     2.2;
activeChannel.gates.beta.q10 =                      2.2;
activeChannel.gates.alpha.equ =                     {'0.0353*(V+27)./(1-exp(-(V+27)/10.2))'};
activeChannel.gates.beta.equ =                      {'0.000883*(-(V+34))./(1-exp((V+34)./10))'};