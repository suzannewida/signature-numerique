% clé publique
p= [173  156   39   16  116  184  166   16   74   61  180  163  192   72   48  108]
 q= [160  142  187  128  144  150  236 235  120  180  134   44   80 189  213   81]
n = p.*q
e =  65537;
% clé privé 
%d = inverse_mod(e,lcm(p-1,q-1));
[w,k,r]= gcd(e,lcm(q-1,p-1));
if w==1  % The inverse of a(mod b) exists only if gcd(a,b)=1
    d = mod(k,lcm(p-1,q-1));
else disp('Modular multiplicative inverse does not exist for these values')
end

%d1 = inverse_mod(e,p-1);
[w,k,r]= gcd(e,p-1);
if w==1  % The inverse of a(mod b) exists only if gcd(a,b)=1
    d1 = mod(k,p-1);
else disp('Modular multiplicative inverse does not exist for these values')
end

%d2 = inverse_mod(e,q-1);
[w,k,r]= gcd(e,q-1);
if w==1  % The inverse of a(mod b) exists only if gcd(a,b)=1
    d2 = mod(k,q-1);
else disp('Modular multiplicative inverse does not exist for these values')
end
%fonction de la signature

m =100
%[r, u, v] = gcd(p,q);
%up = mod(u.*p,q);
%s1 = powermod(m,d1,p);
%s1 = mod(m.^d1,p);
%s2 = powermod(m,d2,q);
%s2 = mod(m.^d2,q);
%iota =  mod(s1 + up.*(s2-s1),n)

signe =  mod(m.^d,(p-1).*(q-1))
design= mod(signe.^e,(p-1).*(q-1))
sha256hasher = System.Security.Cryptography.SHA256Managed;
chiffrement = uint8(sha256hasher.ComputeHash(uint8(signe))) 
