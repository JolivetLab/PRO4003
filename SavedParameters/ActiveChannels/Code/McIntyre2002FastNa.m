function activeChannel = McIntyre2002FastNa()

activeChannel.channames =                           'Fast Na+ (McIntyre 2002)';
activeChannel.cond.value.ref =                      50;
activeChannel.cond.value.vec =                      [];
activeChannel.cond.units =                          {2, 'mS', 'mm', [1, -2]};
activeChannel.erev.value =                          60;
activeChannel.erev.units =                          {1, 'mV', 1};
activeChannel.gates.temp =                          36;
activeChannel.gates.number =                        2;
activeChannel.gates.label =                         {'m', 'h'};
activeChannel.gates.numbereach =                    [3, 1];
activeChannel.gates.alpha.q10 =                     [2.2, 2.9];
activeChannel.gates.beta.q10 =                      [2.2, 2.9];
activeChannel.gates.alpha.equ =                     {'6.57*(V+20.4)./(1-exp(-(V+20.4)./10.3))', ...
                                                         '0.34*(-(V+114))./(1-exp((V+114)./11))'};
activeChannel.gates.beta.equ =                      {'0.304*(-(V+25.7))./(1-exp((V+25.7)./9.16))', ...
                                                         '12.6./(1+exp(-(V+31.8)./13.4))'};