# test_ddh.py
import ddhlib

# Create using new_ddh_precomp
ddh = ddhlib.Ddh.new_ddh_precomp(2, 1024, 100)

# Generate master keys
sk, pk = ddh.generate_master_keys()

print("Secret keys:", sk)
print("Public keys:", pk)
