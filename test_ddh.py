from ddhlib import Ddh
from random import randint

# === Step 1: Setup FE Scheme ===
# Using safe precomputed params (1536-bit for example)
l = 3  # vector length
bound = 10  # absolute value limit for each vector element

fe = Ddh.new_ddh_precomp(l, 1536, bound)

# === Step 2: Generate Keys ===
msk, mpk = Ddh.generate_master_keys(fe)
print("Master Secret Key:", msk)
print("Master Public Key:", mpk)

# === Step 3: Create plaintext vector x and policy vector y ===
x = [randint(-bound, bound) for _ in range(l)]
y = [randint(-bound, bound) for _ in range(l)]

print("Plaintext x:", x)
print("Policy vector y:", y)

# === Step 4: Encrypt x ===
ciphertext = Ddh.encrypt(fe, x, mpk)
print("Ciphertext:", ciphertext)

# === Step 5: Derive functional key for y ===
key = Ddh.derive_key(fe, msk, y)
print("Derived Key (⟨msk, y⟩):", key)

# === Step 6: Decrypt to get ⟨x, y⟩ ===
decrypted_inner_product = Ddh.decrypt(fe, ciphertext, key, y)
print("Decrypted inner product:", decrypted_inner_product)

# === Step 7: Verify correctness ===
expected = sum([x[i] * y[i] for i in range(l)])
print("Expected inner product:", expected)

assert decrypted_inner_product == expected, "❌ FE decryption failed."
print("✅ FE decryption successful.")
