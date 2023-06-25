
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
$\delta \to $N(0.005,0.01)\
$R \to$ IG(5,0.05)\
$\eta_{0} \to$ N([0 0], 0.01*I(2))
