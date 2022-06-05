function activeChannel = McIntyre2002SlowK()

activeChannel.channames =                           'Slow K+ (McIntyre 2002)';
activeChannel.cond.value.ref =                      0.8;
activeChannel.cond.value.vec =                      [];
activeChannel.cond.units =                          {2, 'mS', 'mm', [1, -2]};
activeChannel.erev.value =                          -90;
activeChannel.erev.units =                          {1, 'mV', 1};
activeChannel.gates.temp =                          36;
activeChannel.gates.number =                        1;
activeChannel.gates.label =                         {'s'};
activeChannel.gates.numbereach =                    1;
activeChannel.gates.alpha.q10 =                     3;
activeChannel.gates.beta.q10 =                      3;
activeChannel.gates.alpha.equ =                     {'0.3./(1+exp((V+53)./-5))'};
activeChannel.gates.beta.equ =                      {'0.03./(1+exp((V+90)./-1))'};