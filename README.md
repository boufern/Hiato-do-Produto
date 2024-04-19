
## Hiato do Produto (projeto em desenvolvimento)
-------------------------------------------------------------------------------------
 
Estimação bayesiana de um modelo de espaço-estado capaz de estimar o hiato de produto e
o produto potencial da economia brasileira.

### Modelo dinâmico com parâmetros constantes:


$y_{t} = \mu_{t} + \eta_{t}$\
$\mu_{t} = \delta + \mu_{t-1} + w_{t}$                          , $w(t) \to N(0,Q)$ iid\
$\eta_{t} = \phi_{1}\eta_{t-1} + \phi_{2}\eta_{t-2} + v_{t}$    , $v(t) \to N(0,R)$ iid

onde:\
$y_{t}$ -- log do produto\
$\mu(t)$ -- produto potencial\
$\eta(t)$ -- hiato do produto
 
### Priors e Valores Iniciais

$\phi_{1} \to$ N(0.7,0.1)\
$\phi_{2} \to$ N(0,0.05)\
$\delta \to$ N(0.005,0.01)\
$Q \to$ IG(5,0.05)\
$R \to$ IG(5,0.05)\
$\eta_{0} \to$ N([0 0], 0.01 $\times$ I(2))\

As priors refletem que é provável que $\eta_{t}$ seja estacionário, que o crescimento do produto seja mais ou menos 2% ao ano com incerteza substancial em torno de sua média. As priors das variâncias são definidas com razoável incerteza de modo que suas médias sejam a metade da variância total do log do produto. Os estado inicial $\eta_{0}$ = [ $\eta_{t}$ $\eta_{t-1}$ ]' impõe a crença de que a economia estava próxima à tendência no começo da amostra de dados.


