# Deploying SeqLit Online (Free Hosting & Free Domain)

This guide walks you through publishing **SeqLit** on the web with a **free public URL / domain** and **free cloud hosting** directly connected to your GitHub repository (`https://github.com/rizwanahmed2603-byte/SQLit.git`).

Whenever you push new changes to `main` on GitHub, your live website will **automatically rebuild and redeploy**.

---

## Option 1: Render (Recommended - 2 Minutes)
**Free Domain:** `https://seqlit.onrender.com` (or your chosen name)  
**Cost:** 100% Free Tier (Free Web Service)  
**Auto-deploy:** Yes (connected to your GitHub repository)

### Steps:
1. Go to **[render.com](https://render.com)** and sign up or log in with your **GitHub account**.
2. Click the **"New +"** button in the dashboard and select **"Web Service"**.
3. Under *Connect a repository*, find and select your repo: **`SQLit`** (or paste `https://github.com/rizwanahmed2603-byte/SQLit.git`).
4. Configure the service settings (Render auto-detects most settings from our repository configuration):
   - **Name:** `seqlit` (your site will be `https://seqlit.onrender.com`)
   - **Region:** Choose closest to your users (e.g., Singapore, Frankfurt, Oregon, Ohio)
   - **Branch:** `main`
   - **Root Directory:** *(leave blank)*
   - **Runtime:** `Python 3`
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `gunicorn app:app --timeout 120`
   - **Instance Type:** Select **Free**
5. Click **"Deploy Web Service"** at the bottom.
6. Wait 2–3 minutes for the build to finish. Once done, Render gives you your free live URL:
   ```
   https://seqlit.onrender.com
   ```

> [!TIP]
> The free tier on Render automatically spins down after 15 minutes of inactivity. When a new visitor opens your link, it takes ~30 seconds to wake up.

---

## Option 2: Hugging Face Spaces (Ideal for Science & AI Tools)
**Free Domain:** `https://<your-username>-seqlit.hf.space`  
**Cost:** 100% Free Forever (No credit card required)  
**Specs:** 2 vCPU, 16 GB RAM (generous for heavy alignments)

### Steps:
1. Go to **[huggingface.co](https://huggingface.co)** and log in or create a free account.
2. In the top-right profile menu, click **"New Space"**.
3. Configure your Space:
   - **Space Name:** `seqlit`
   - **License:** `mit` or `apache-2.0`
   - **Space SDK:** Select **Docker** -> **Blank**
   - **Space Hardware:** Free (2 vCPU · 16 GB RAM)
   - **Visibility:** Public
4. Click **"Create Space"**.
5. To link your GitHub repo:
   - In your newly created Space, go to **Settings** -> **Repository settings**.
   - You can push your GitHub code directly to the Space's Git remote or set up GitHub Actions to auto-sync.
   - Alternatively, use the Hugging Face web UI to upload files or clone from GitHub.
6. Hugging Face automatically detects our [`Dockerfile`](Dockerfile) and starts SeqLit on port 7860!

---

## Option 3: Koyeb (Always-On Free Nano Tier)
**Free Domain:** `https://seqlit.koyeb.app`  
**Cost:** 100% Free "Eco" tier

### Steps:
1. Go to **[koyeb.com](https://koyeb.com)** and log in with GitHub.
2. Click **"Create App"**.
3. Select **GitHub** as the deployment method and choose `SQLit`.
4. Koyeb automatically reads our `Dockerfile` or `Procfile`.
5. Click **"Deploy"**. Your app will be live at `https://<your-app>.koyeb.app`.

---

## Adding a Custom Domain (Optional)
If you later register a custom domain (e.g., `seqlit.bio` or a free subdomain like `seqlit.duckdns.org` or `seqlit.is-a.dev`):
1. In your **Render** dashboard, open your service and click **Settings** -> **Custom Domains**.
2. Click **Add Custom Domain** and enter your domain name.
3. Add the provided CNAME record in your domain's DNS manager.
4. Render automatically provisions a free Let's Encrypt SSL certificate for HTTPS!
